//! Python integer true division for the DAS TSV's unbounded coordinates.
//! Integer quotient/remainder arithmetic rounds once to binary64, including
//! ties, subnormal values and signed underflow. The existing bounded metadata
//! and count helper remains a separate, unchanged implementation.
use num_bigint::BigInt;

#[derive(Debug, PartialEq)]
pub(crate) enum Failure {
    ZeroDivision,
    Overflow,
}
impl Failure {
    pub(crate) fn category(&self) -> &'static str {
        match self {
            Self::ZeroDivision => "ZeroDivisionError",
            Self::Overflow => "OverflowError",
        }
    }
}

pub(crate) fn divide(numerator: &BigInt, denominator: &BigInt) -> Result<f64, Failure> {
    use num_bigint::Sign;
    let sign = u64::from((numerator.sign() == Sign::Minus) != (denominator.sign() == Sign::Minus));
    let a = numerator.magnitude();
    let b = denominator.magnitude();
    if b.bits() == 0 {
        return Err(Failure::ZeroDivision);
    }
    let signed_zero = f64::from_bits(sign << 63);
    if a.bits() == 0 {
        return Ok(signed_zero);
    }
    let mut exponent = i128::from(a.bits()) - i128::from(b.bits());
    // A bit-length difference is at most one above floor(log2(a / b)).
    // These guards avoid allocating a shift proportional to a huge exponent.
    if exponent > 1024 {
        return Err(Failure::Overflow);
    }
    if exponent < -1075 {
        return Ok(signed_zero);
    }
    let below = if exponent >= 0 {
        a < &(b << exponent as usize)
    } else {
        &(a << (-exponent) as usize) < b
    };
    if below {
        exponent -= 1;
    }
    if exponent > 1023 {
        return Err(Failure::Overflow);
    }
    if exponent < -1075 {
        return Ok(signed_zero);
    }
    let subnormal = exponent < -1022;
    let shift = if subnormal { 1074 } else { 52 - exponent };
    let (scaled_a, scaled_b) = if shift >= 0 {
        (a << shift as usize, b.clone())
    } else {
        (a.clone(), b << (-shift) as usize)
    };
    let quotient = &scaled_a / &scaled_b;
    let remainder = &scaled_a % &scaled_b;
    // The scaling above bounds the quotient to 53 bits, regardless of either
    // input's decimal size. Subnormal zero has no nonzero limbs.
    let digits = quotient.to_u64_digits();
    debug_assert!(digits.len() <= 1);
    let mut significand = digits.first().copied().unwrap_or(0);
    let twice_remainder = remainder << 1usize;
    if twice_remainder > scaled_b || (twice_remainder == scaled_b && significand & 1 != 0) {
        significand += 1;
    }
    if subnormal {
        // Carry to 2^52 is precisely the smallest normal value's bit pattern.
        return Ok(f64::from_bits((sign << 63) | significand));
    }
    if significand == 1 << 53 {
        significand >>= 1;
        exponent += 1;
        if exponent > 1023 {
            return Err(Failure::Overflow);
        }
    }
    Ok(f64::from_bits(
        (sign << 63) | (((exponent + 1023) as u64) << 52) | (significand & ((1 << 52) - 1)),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary64_ties_and_subnormal_boundaries() {
        let one = BigInt::from(1);
        let denominator = &one << 1075usize;
        assert_eq!(divide(&one, &denominator).unwrap().to_bits(), 0);
        assert_eq!(divide(&BigInt::from(3), &denominator).unwrap().to_bits(), 2);
        assert_eq!(divide(&(-&one), &denominator).unwrap().to_bits(), 1 << 63);
        assert_eq!(divide(&one, &(&denominator - 1)).unwrap().to_bits(), 1);
        assert_eq!(divide(&one, &(&denominator + 1)).unwrap().to_bits(), 0);
        let scale = &one << 53usize;
        assert_eq!(divide(&(&scale + 1), &scale).unwrap(), 1.0);
        assert_eq!(
            divide(&(&scale + 3), &scale).unwrap().to_bits(),
            1.0f64.to_bits() + 2
        );
    }

    #[test]
    fn signed_zero_zero_division_and_overflow() {
        let zero = BigInt::from(0);
        let one = BigInt::from(1);
        assert_eq!(divide(&zero, &(-&one)).unwrap().to_bits(), 1 << 63);
        assert_eq!(divide(&one, &zero), Err(Failure::ZeroDivision));
        assert_eq!(divide(&(&one << 1024usize), &one), Err(Failure::Overflow));
        let boundary = ((&one << 54usize) - 1) << 970usize;
        assert_eq!(divide(&boundary, &one), Err(Failure::Overflow));
        assert_eq!(divide(&(&boundary - 1), &one).unwrap(), f64::MAX);
    }
}
