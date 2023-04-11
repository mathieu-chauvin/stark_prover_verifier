// implement polynomials using the FieldElement type

use crate::field::FieldElement;

// Poly is a polynomial of degree n
#[derive(Copy, Clone)]
pub struct Poly {
    coeffs: [FieldElement; 256],
    // type can be coeffs or lagrange
}



// implement add, mul
// implement eval
// implement div

// implement the Lagrange interpolation algorithm
// https://en.wikipedia.org/wiki/Lagrange_polynomial



// implement the FFT algorithm
// https://en.wikipedia.org/wiki/Fast_Fourier_transform

impl Poly {
    
    fn new(coeffs: [FieldElement; 256]) -> Poly {
        Poly { coeffs: coeffs }
    }

    fn add(&self, other: Poly) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            result[i] = self.coeffs[i] + other.coeffs[i];
        }
        Poly::new(result)
    }

    fn sub(&self, other: Poly) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            result[i] = self.coeffs[i]- other.coeffs[i];
        }
        Poly::new(result)
    }

    fn eval(&self, x: FieldElement) -> FieldElement {
        let mut result = FieldElement::new(0);
        for i in 0..256 {
            result = result+(self.coeffs[i]*x.pow(i as u64));
        }
        result
    }

    //interpolation
    
}

//tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let a = Poly::new([FieldElement::new(1);256]);
        let b = Poly::new([FieldElement::new(2); 256]);
        let c = a.add(b);
        assert_eq!(c.coeffs[0], FieldElement::new(3));
    }

    #[test]
    fn test_sub() {
        let a = Poly::new([FieldElement::new(1);256]);
        let b = Poly::new([FieldElement::new(2); 256]);
        let c = a.sub(b);
        assert_eq!(c.coeffs[0], FieldElement::new(0xffffffff00000000));
    }

    #[test]
    fn test_eval() {
        let mut a = Poly::new([FieldElement::new(0);256]);
        a.coeffs[0] = FieldElement::new(2);
        a.coeffs[1] = FieldElement::new(1);

        assert_eq!(a.eval(FieldElement::new(1)), FieldElement::new(3));
    }
}