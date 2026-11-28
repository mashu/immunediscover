//! Telling from the header of a genome record which chromosome it belongs
//! to.
//!
//! Genome files name their records differently. UCSC writes `chr14` and
//! `chr14_KI270846v1_alt`; NCBI writes an accession followed by
//! `Homo sapiens chromosome 14, GRCh38.p14 Primary Assembly`; Ensembl
//! writes `14 dna:chromosome ...`. All three are understood here.

/// The chromosome a record belongs to (`"14"`, `"X"`), or `None` for
/// records that belong to none (unplaced scaffolds, the mitochondrion of
/// an NCBI file).
pub fn chromosome_of(header: &str) -> Option<String> {
    let name = header.split_ascii_whitespace().next().unwrap_or("");
    let label = |text: &str| {
        let end = text
            .find(|c: char| !c.is_ascii_alphanumeric())
            .unwrap_or(text.len());
        (end > 0).then(|| text[..end].to_ascii_uppercase())
    };
    // UCSC: chr14, chr14_GL000009v2_random, chrUn_KI270302v1.
    if name.len() > 3 && name[..3].eq_ignore_ascii_case("chr") {
        return label(&name[3..]).filter(|found| found != "UN");
    }
    // Ensembl: the name is the chromosome itself.
    let plain = name.to_ascii_uppercase();
    let is_number =
        plain.len() <= 2 && !plain.is_empty() && plain.bytes().all(|b| b.is_ascii_digit());
    if is_number || matches!(plain.as_str(), "X" | "Y" | "MT") {
        return Some(plain);
    }
    // NCBI: "... chromosome 14, ..." in the description.
    let lower = header.to_ascii_lowercase();
    let at = lower.find("chromosome ")? + "chromosome ".len();
    label(&header[at..])
}

/// The chromosome a user asks for, from `14`, `chr14` or `chromosome 14`.
pub fn wanted_chromosome(text: &str) -> String {
    let trimmed = text.trim();
    let lower = trimmed.to_ascii_lowercase();
    let rest = ["chromosome", "chr"]
        .iter()
        .find_map(|prefix| lower.strip_prefix(prefix).map(|_| &trimmed[prefix.len()..]))
        .unwrap_or(trimmed);
    rest.trim().to_ascii_uppercase()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_three_common_ways_of_naming_records_are_understood() {
        let cases = [
            ("chr14", Some("14")),
            ("chr14_KI270846v1_alt", Some("14")),
            ("chr14_GL000009v2_random", Some("14")),
            ("chrX", Some("X")),
            ("chrUn_KI270302v1", None),
            ("chrM", Some("M")),
            (
                "NC_000014.9 Homo sapiens chromosome 14, GRCh38.p14 Primary Assembly",
                Some("14"),
            ),
            (
                "NT_187600.1 Homo sapiens chromosome 14 genomic contig, GRCh38.p14 alternate locus group ALT_REF_LOCI_1 HSCHR14_3_CTG1",
                Some("14"),
            ),
            (
                "NC_060938.1 Homo sapiens isolate CHM13 chromosome 14, alternate assembly T2T-CHM13v2.0",
                Some("14"),
            ),
            ("NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly", Some("1")),
            ("NT_187361.1 Homo sapiens unplaced genomic scaffold, GRCh38.p14 Primary Assembly HSCHRUN_RANDOM_CTG1", None),
            ("NC_012920.1 Homo sapiens mitochondrion, complete genome", None),
            ("14 dna:chromosome chromosome:GRCh38:14:1:107043718:1 REF", Some("14")),
            ("X dna:chromosome chromosome:GRCh38:X:1:156040895:1 REF", Some("X")),
            ("KI270728.1 dna:scaffold scaffold:GRCh38:KI270728.1:1:1872759:1 REF", None),
        ];
        for (header, expected) in cases {
            assert_eq!(chromosome_of(header).as_deref(), expected, "{header}");
        }
    }

    #[test]
    fn chromosome_one_is_not_chromosome_fourteen() {
        let fourteen = "NC_000014.9 Homo sapiens chromosome 14, GRCh38.p14 Primary Assembly";
        assert_ne!(chromosome_of(fourteen).unwrap(), wanted_chromosome("1"));
        assert_eq!(
            chromosome_of(fourteen).unwrap(),
            wanted_chromosome("chromosome 14")
        );
    }

    #[test]
    fn the_wanted_chromosome_can_be_written_in_several_ways() {
        for text in ["14", "chr14", "Chr14", "chromosome 14", " 14 "] {
            assert_eq!(wanted_chromosome(text), "14", "{text}");
        }
        assert_eq!(wanted_chromosome("x"), "X");
        assert_eq!(wanted_chromosome("chromosome 22"), "22");
    }
}
