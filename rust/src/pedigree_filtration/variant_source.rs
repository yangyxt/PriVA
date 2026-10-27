//! Source rules do not repeat the dispatcher's proband gate on direct calls.
use super::SourceResult;
use super::genotype::Genotype;

pub fn chromosome_type(chromosome: &str) -> &'static str {
    match chromosome.replace("chr", "").to_uppercase().as_str() {
        "X" => "X",
        "Y" => "Y",
        _ => "AUTOSOMAL",
    }
}

pub fn autosomal(
    proband: &Genotype,
    father: &Genotype,
    mother: &Genotype,
) -> SourceResult<&'static str> {
    if !father.exists()? || !mother.exists()? {
        return Ok(if father.has_variant()? || mother.has_variant()? {
            "INHERITED"
        } else {
            "UNKNOWN"
        });
    }
    // The original evaluates this unused local AFTER both exists checks.
    // Preserve its possible error/evaluation order instead of optimizing it out.
    let _proband_has_variant = proband.has_variant()?;
    let paternal = father.has_variant()?;
    let maternal = mother.has_variant()?;
    Ok(if !paternal && !maternal {
        "DENOVO"
    } else if paternal && !maternal {
        "PATERNAL"
    } else if !paternal && maternal {
        "MATERNAL"
    } else if paternal && maternal {
        "INHERITED"
    } else {
        "UNKNOWN"
    })
}

pub fn chr_x(
    proband: &Genotype,
    sex: &str,
    father: &Genotype,
    mother: &Genotype,
) -> SourceResult<&'static str> {
    match sex {
        "1" => Ok(if !mother.exists()? {
            "UNKNOWN"
        } else if mother.has_variant()? {
            "MATERNAL"
        } else {
            "DENOVO"
        }),
        "2" => autosomal(proband, father, mother),
        _ => Ok("UNKNOWN"),
    }
}

pub fn chr_y(sex: &str, father: &Genotype) -> SourceResult<&'static str> {
    Ok(if sex != "1" || !father.exists()? {
        "UNKNOWN"
    } else if father.has_variant()? {
        "PATERNAL"
    } else {
        "DENOVO"
    })
}

pub fn source(
    proband: &Genotype,
    chromosome_type: &str,
    sex: &str,
    father: &Genotype,
    mother: &Genotype,
) -> SourceResult<&'static str> {
    match chromosome_type {
        "AUTOSOMAL" => autosomal(proband, father, mother),
        "X" => chr_x(proband, sex, father, mother),
        "Y" => chr_y(sex, father),
        _ => Ok("UNKNOWN"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pedigree_filtration::genotype::Allele;

    #[test]
    fn unused_proband_check_remains_after_parent_completeness() {
        let invalid = Genotype::Scalar("int");
        let reference = Genotype::Sequence(vec![Allele::Integer("0".into())]);
        assert!(autosomal(&invalid, &reference, &reference).is_err());
        assert_eq!(
            autosomal(&invalid, &Genotype::None, &reference),
            Ok("UNKNOWN")
        );
        // Direct calls deliberately do not reapply the proband ALT1 gate.
        assert_eq!(
            autosomal(&Genotype::None, &reference, &reference),
            Ok("DENOVO")
        );
    }
}
