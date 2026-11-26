//! Fitting the D-gene model and scanning reads with it.

use super::model::{
    encode, DurationModel, EmpiricalDuration, IidEmission, Pwm, TruncatedNormalDuration,
};
use crate::error::Result;

/// Additive smoothing used for every model component.
pub const PSEUDOCOUNT: f64 = 0.01;

const NONAMER: usize = 9;
const SPACER: usize = 12;
const HEPTAMER: usize = 7;
/// Length of one complete recombination signal.
const SIGNAL: usize = NONAMER + SPACER + HEPTAMER;

/// One training example: a known D allele with its two signals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrainingExample {
    pub pre_nonamer: String,
    pub pre_spacer: String,
    pub pre_heptamer: String,
    pub gene: String,
    pub post_heptamer: String,
    pub post_spacer: String,
    pub post_nonamer: String,
}

/// The D gene with flanking recombination signals.
#[derive(Debug, Clone, PartialEq)]
pub struct DGeneRssHsmm<D: DurationModel = EmpiricalDuration> {
    pub pre_nonamer: Pwm<NONAMER>,
    pub pre_spacer: Pwm<SPACER>,
    pub pre_heptamer: Pwm<HEPTAMER>,
    pub gene: IidEmission,
    pub duration: D,
    pub post_heptamer: Pwm<HEPTAMER>,
    pub post_spacer: Pwm<SPACER>,
    pub post_nonamer: Pwm<NONAMER>,
}

impl<D: DurationModel> DGeneRssHsmm<D> {
    /// Fit the emission models from training examples and use `duration`
    /// for the gene length.
    pub fn fit_with(examples: &[TrainingExample], duration: D) -> Result<Self> {
        let column = |pick: fn(&TrainingExample) -> &String| -> Vec<&str> {
            examples.iter().map(|e| pick(e).as_str()).collect()
        };
        Ok(DGeneRssHsmm {
            pre_nonamer: Pwm::fit(&column(|e| &e.pre_nonamer), PSEUDOCOUNT)?,
            pre_spacer: Pwm::fit(&column(|e| &e.pre_spacer), PSEUDOCOUNT)?,
            pre_heptamer: Pwm::fit(&column(|e| &e.pre_heptamer), PSEUDOCOUNT)?,
            gene: IidEmission::fit(&column(|e| &e.gene), PSEUDOCOUNT),
            duration,
            post_heptamer: Pwm::fit(&column(|e| &e.post_heptamer), PSEUDOCOUNT)?,
            post_spacer: Pwm::fit(&column(|e| &e.post_spacer), PSEUDOCOUNT)?,
            post_nonamer: Pwm::fit(&column(|e| &e.post_nonamer), PSEUDOCOUNT)?,
        })
    }
}

impl DGeneRssHsmm<EmpiricalDuration> {
    /// Fit all components from training examples; gene lengths outside
    /// `min_gene..=max_gene` do not inform the duration model.
    pub fn fit(examples: &[TrainingExample], min_gene: usize, max_gene: usize) -> Result<Self> {
        let genes: Vec<&str> = examples.iter().map(|e| e.gene.as_str()).collect();
        let duration = EmpiricalDuration::fit(&genes, min_gene, max_gene, PSEUDOCOUNT.max(1e-6))?;
        Self::fit_with(examples, duration)
    }
}

impl DGeneRssHsmm<TruncatedNormalDuration> {
    /// Fit with a discretised truncated normal gene-length model (the
    /// original's `duration = :normal`).
    pub fn fit_normal(
        examples: &[TrainingExample],
        min_gene: usize,
        max_gene: usize,
        mean: f64,
        sd: f64,
    ) -> Result<Self> {
        Self::fit_with(
            examples,
            TruncatedNormalDuration::new(min_gene, max_gene, mean, sd)?,
        )
    }
}

/// The most probable placement of the model on a read.
///
/// All positions are 1-based and inclusive.
#[derive(Debug, Clone, PartialEq)]
pub struct Detection {
    pub prefix_start: usize,
    pub prefix_end: usize,
    pub gene_start: usize,
    pub gene_end: usize,
    pub suffix_start: usize,
    pub suffix_end: usize,
    /// Log-probability of the best placement.
    pub log_path_prob: f64,
    /// Log of the summed probability of all placements.
    pub log_total_prob: f64,
    /// Share of the total probability carried by the best placement.
    pub posterior_prob: f64,
}

/// `log(exp(a) + exp(b))` without overflow.
#[inline]
fn log_add_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let max = a.max(b);
    max + (a.min(b) - max).exp().ln_1p()
}

impl<D: DurationModel> DGeneRssHsmm<D> {
    /// Log-probability of the gene state emitting `duration` bases from `start` (0-based).
    #[inline]
    fn score_gene(&self, observations: &[Option<u8>], start: usize, duration: usize) -> f64 {
        let mut total = self.duration.log_pmf(duration);
        for observation in &observations[start..start + duration] {
            match observation {
                Some(base) => total += self.gene.log_p[*base as usize],
                None => return f64::NEG_INFINITY,
            }
        }
        total
    }

    /// Evaluate every placement of the model on `sequence`; return the best
    /// one together with the total probability, or `None` if nothing fits.
    pub fn scan(&self, sequence: &str) -> Option<Detection> {
        let observations = encode(sequence);
        let len = observations.len();
        let (min_gene, max_gene) = (self.duration.min_len(), self.duration.max_len());
        let shortest = SIGNAL + min_gene + SIGNAL;
        if len < shortest {
            return None;
        }

        let mut best_log = f64::NEG_INFINITY;
        let mut total_log = f64::NEG_INFINITY;
        // (signal start, gene start, gene length), 0-based.
        let mut best: Option<(usize, usize, usize)> = None;

        for signal_start in 0..=(len - shortest) {
            let mut upstream = self.pre_nonamer.score(&observations, signal_start);
            upstream += self.pre_spacer.score(&observations, signal_start + NONAMER);
            upstream += self
                .pre_heptamer
                .score(&observations, signal_start + NONAMER + SPACER);
            if upstream == f64::NEG_INFINITY {
                continue;
            }
            let gene_start = signal_start + SIGNAL;
            for duration in min_gene..=max_gene {
                let downstream_start = gene_start + duration;
                if downstream_start + SIGNAL > len {
                    break;
                }
                let gene = self.score_gene(&observations, gene_start, duration);
                if gene == f64::NEG_INFINITY {
                    continue;
                }
                let downstream = self.post_heptamer.score(&observations, downstream_start)
                    + self
                        .post_spacer
                        .score(&observations, downstream_start + HEPTAMER)
                    + self
                        .post_nonamer
                        .score(&observations, downstream_start + HEPTAMER + SPACER);
                if downstream == f64::NEG_INFINITY {
                    continue;
                }
                let total = upstream + gene + downstream;
                total_log = log_add_exp(total_log, total);
                if total > best_log {
                    best_log = total;
                    best = Some((signal_start, gene_start, duration));
                }
            }
        }

        let (signal_start, gene_start, duration) = best?;
        let log_posterior = best_log - total_log;
        Some(Detection {
            prefix_start: signal_start + 1,
            prefix_end: signal_start + SIGNAL,
            gene_start: gene_start + 1,
            gene_end: gene_start + duration,
            suffix_start: gene_start + duration + 1,
            suffix_end: gene_start + duration + SIGNAL,
            log_path_prob: best_log,
            log_total_prob: total_log,
            posterior_prob: if log_posterior.is_finite() {
                log_posterior.exp()
            } else {
                0.0
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PRE_NONAMER: &str = "GGTTTTTGT";
    const PRE_SPACER: &str = "ACTGATCGATCA";
    const PRE_HEPTAMER: &str = "CACTGTG";
    const POST_HEPTAMER: &str = "CACAGTG";
    const POST_SPACER: &str = "TTGACCGATGAC";
    const POST_NONAMER: &str = "ACAAAAACC";

    fn example(gene: &str) -> TrainingExample {
        TrainingExample {
            pre_nonamer: PRE_NONAMER.into(),
            pre_spacer: PRE_SPACER.into(),
            pre_heptamer: PRE_HEPTAMER.into(),
            gene: gene.into(),
            post_heptamer: POST_HEPTAMER.into(),
            post_spacer: POST_SPACER.into(),
            post_nonamer: POST_NONAMER.into(),
        }
    }

    fn read(gene: &str) -> String {
        format!(
            "TTTTTTCCCCCCGGGGGG{PRE_NONAMER}{PRE_SPACER}{PRE_HEPTAMER}{gene}{POST_HEPTAMER}{POST_SPACER}{POST_NONAMER}AAACCCGGGTTTAAA"
        )
    }

    #[test]
    fn locates_a_d_gene_between_its_signals() {
        let genes = [
            "GGTATAGTGGGAGCTAC",
            "AGGATATTGTAGTGGTGGTAGC",
            "GTATTACTATGGTTCGGGGAGTTATTATAAC",
        ];
        let examples: Vec<TrainingExample> = genes.iter().map(|g| example(g)).collect();
        let model = DGeneRssHsmm::fit(&examples, 10, 40).unwrap();
        let novel = "GGTATAGTGGGAGGTAC";
        let sequence = read(novel);
        let detection = model.scan(&sequence).unwrap();
        assert_eq!(
            &sequence[detection.gene_start - 1..detection.gene_end],
            novel
        );
        assert_eq!(
            &sequence[detection.prefix_start - 1..detection.prefix_end].len(),
            &28
        );
        assert_eq!(detection.suffix_end - detection.suffix_start + 1, 28);
        assert!(detection.posterior_prob > 0.9);
        assert!(detection.log_total_prob >= detection.log_path_prob);
    }

    #[test]
    fn short_or_ambiguous_reads_yield_nothing() {
        let model = DGeneRssHsmm::fit(&[example("GGTATAGTGGGAGCTAC")], 10, 40).unwrap();
        assert!(model.scan("ACGT").is_none());
        assert!(model.scan(&"N".repeat(200)).is_none());
    }

    #[test]
    fn log_add_exp_handles_infinities() {
        assert_eq!(log_add_exp(f64::NEG_INFINITY, -3.0), -3.0);
        assert_eq!(log_add_exp(-3.0, f64::NEG_INFINITY), -3.0);
        assert!((log_add_exp(0.0, 0.0) - 2f64.ln()).abs() < 1e-15);
    }
}
