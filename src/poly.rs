// implement polynomials using the FieldElement type

use crate::field::FieldElement;
use crate::field::P;

// implement a polynomial type
// enum 
// 1. coeffs
// 2. lagrange



// Poly is a polynomial of degree n
#[derive(Copy, Clone)]
pub struct Poly {
    coeffs: [FieldElement; 256],
    // type can be coeffs or lagrange
    lagrange : bool

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
        Poly { coeffs: coeffs, lagrange: false }
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

    fn mul(&self, other: FieldElement) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            result[i] = self.coeffs[i]*other;
            
        }
        Poly::new(result)
    }

    fn mul_poly(&self, other: Poly) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            for j in 0..256 {
                if i+j < 256 {
                    result[i+j] = result[i+j]+(self.coeffs[i]*other.coeffs[j]);
                }
            }
            
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

    fn lagrange_interpolation(x: &[FieldElement], y:&[FieldElement]) -> Poly {
        assert_eq!(x.len(), y.len());

        let coeffs = [FieldElement::new(0); 256];
        let mut result = Poly::new(coeffs);
        for i in 0..x.len() {
            let mut coeffs = [FieldElement::new(0); 256];
            coeffs[0] = FieldElement::new(1);
            let mut numerator = Poly::new(coeffs);

            let mut denominator = FieldElement::new(1);
            for j in 0..y.len() {
                if i != j {

                    //X - Xj
                    let mut coeffs = [FieldElement::new(0); 256];
                    coeffs[0] = FieldElement::new(P-(j as u64));
                    coeffs[1] = FieldElement::new(1);
                    let poly = Poly::new(coeffs);
                    
                    numerator = numerator.mul_poly(poly);
                    denominator = denominator*(FieldElement::new(i as u64)-FieldElement::new(j as u64));
                }
            }
            result = result.add(numerator.mul(y[i]).mul(denominator.inv()));
        }
        result
    }


    // fast fourier transform
    // https://en.wikipedia.org/wiki/Fast_Fourier_transform
    

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

    #[test]
    fn test_mul() {
        let a = Poly::new([FieldElement::new(1);256]);
        let b = a.mul(FieldElement::new(2));
        assert_eq!(b.coeffs[0], FieldElement::new(2));
    }

    #[test]
    fn test_mul_poly() {
        let mut coeffs = [FieldElement::new(0); 256];
        coeffs[2] = FieldElement::new(1);
        let a = Poly::new(coeffs);
        let b = Poly::new([FieldElement::new(2);256]);
        let c = a.mul_poly(b);
        assert_eq!(c.coeffs[0], FieldElement::new(0));
        assert_eq!(c.coeffs[2], FieldElement::new(2));
    }

    #[test]
    fn test_lagrange_interpolation() {
        let x = [FieldElement::new(0), FieldElement::new(1), FieldElement::new(2)];
        let y = [FieldElement::new(2), FieldElement::new(3), FieldElement::new(6)];

        //let expected_coeffs = [FieldElement::new(2), FieldElement::new(0), FieldElement::new(1)];
        let mut expected_coeffs = [FieldElement::new(0); 256];
        expected_coeffs[0] = FieldElement::new(2);
        expected_coeffs[1] = FieldElement::new(0);
        expected_coeffs[2] = FieldElement::new(1);


        let p = Poly::lagrange_interpolation(&x, &y);
        assert_eq!(p.coeffs, expected_coeffs);

    }
}