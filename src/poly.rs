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
    pub coeffs: [FieldElement; 256],
    // type can be coeffs or lagrange
    pub lagrange : bool

}


impl Poly {
    
    pub fn new(coeffs: [FieldElement; 256]) -> Poly {
        Poly { coeffs: coeffs, lagrange: false }
    }

    pub fn add(&self, other: Poly) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            result[i] = self.coeffs[i] + other.coeffs[i];
        }
        Poly::new(result)
    }

    pub fn sub(&self, other: Poly) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            result[i] = self.coeffs[i]- other.coeffs[i];
        }
        Poly::new(result)
    }

    pub fn mul(&self, other: FieldElement) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            result[i] = self.coeffs[i]*other;
            
        }
        Poly::new(result)
    }

    pub fn mul_poly(&self, other: Poly) -> Poly {
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


    pub fn eval(&self, x: FieldElement) -> FieldElement {
        let mut result = FieldElement::new(0);
        for i in 0..256 {
            result = result+(self.coeffs[i]*x.pow(i as u64));
        }
        result
    }

    pub fn lagrange_interpolation(x: &[FieldElement], y:&[FieldElement]) -> Poly {
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
    
    pub fn fft(vals: &[FieldElement], root_of_unity: &FieldElement) -> Vec<FieldElement> {
        let n = vals.len();
        if n == 1 {
            return vals.to_vec();
        }

        let root_of_unity_sqr = (*root_of_unity)*(*root_of_unity);
        let l = Self::fft(&vals.iter().step_by(2).map(|&x| x).collect::<Vec<FieldElement>>(), &root_of_unity_sqr);
        let r = Self::fft(&vals.iter().skip(1).step_by(2).map(|&x| x).collect::<Vec<FieldElement>>(), &root_of_unity_sqr);

        let mut o = vec![FieldElement::new(0); n];
        let mut power_of_root_of_unity = FieldElement::new(1);
            
        for i in 0..n/2 {
            let y_times_root = r[i]*power_of_root_of_unity;
            o[i] = l[i]+(y_times_root);
            o[i+n/2] = l[i]-(y_times_root);
            power_of_root_of_unity = power_of_root_of_unity*(*root_of_unity);
        }
        o
    }

    // Naive FFT function, for test purposes
    fn naive_fft(vals: &[FieldElement], root_of_unity: &FieldElement) -> Vec<FieldElement> {
        let n = vals.len();
        let mut o = vec![FieldElement::new(0); n];
        let mut power_of_root_of_unity = FieldElement::new(1);

        for i in 0..n {
            for j in 0..n {
                o[i] = o[i]+(vals[j]*power_of_root_of_unity.pow(j as u64));
            }
            power_of_root_of_unity = power_of_root_of_unity*(*root_of_unity);
        }
        o
    }

    // Inverse FFT function
    pub fn inv_fft(vals: &[FieldElement], root_of_unity: &FieldElement) -> Vec<FieldElement> {
        let n = vals.len();

        // Inverse FFT
        let invlen = FieldElement::new(vals.len() as u64).inv();
        let inv_root_of_unity = root_of_unity.inv();
        print!("invlen: {:?}" , invlen);
        let poly = Self::fft(&vals, &inv_root_of_unity);
        print!("poly: {:?}" , poly);
        // for each element of poly, multiply by invlen
        let inv_poly: Vec<FieldElement> = poly.iter().map(|x| (*x * invlen)).collect();
        inv_poly

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

    #[test]
    fn test_fft() {
        let mut coeffs = [FieldElement::new(0);8];
        coeffs[1] = FieldElement::new(1);

        let n = FieldElement::nth_root_of_unity(8);
        let fft =Poly::fft(&coeffs, &n);
        let mut arr = [FieldElement::new(0); 256];
        let len = fft.len();
        arr[..len].copy_from_slice(&fft); // copy the elements from the vector to the array
        
        let expected_coeffs = [FieldElement::new(1),n, n.pow(2), n.pow(3), n.pow(4),n.pow(5),n.pow(6),n.pow(7)];
        
        assert_eq!(arr[..len], expected_coeffs);
    }

    #[test]
    fn test_naive_fft() {
        let mut coeffs = [FieldElement::new(0);8];
        coeffs[1] = FieldElement::new(1);

        let n = FieldElement::nth_root_of_unity(8);
        let fft =Poly::naive_fft(&coeffs, &n);
        let mut arr = [FieldElement::new(0); 256];
        let len = fft.len();
        arr[..len].copy_from_slice(&fft); // copy the elements from the vector to the array
        
        let expected_coeffs = [FieldElement::new(1),n, n.pow(2), n.pow(3), n.pow(4),n.pow(5),n.pow(6),n.pow(7)];
        
        assert_eq!(arr[..len], expected_coeffs);
    }


    #[test]
    fn test_fft_2() {
        let mut coeffs = [FieldElement::new(0);8];
        coeffs[1] = FieldElement::new(1);
        coeffs[0] = FieldElement::new(P-1);

        let n = FieldElement::nth_root_of_unity(8);
        let fft =Poly::fft(&coeffs, &n);
        let mut arr = [FieldElement::new(0); 256];
        let len = fft.len();
        arr[..len].copy_from_slice(&fft); // copy the elements from the vector to the array
        
        let expected_coeffs = [FieldElement::new(0),n-FieldElement::new(1), n.pow(2)-FieldElement::new(1), n.pow(3)-FieldElement::new(1), n.pow(4)-FieldElement::new(1),n.pow(5)-FieldElement::new(1),n.pow(6)-FieldElement::new(1),n.pow(7)-FieldElement::new(1)];
        
        assert_eq!(arr[..len], expected_coeffs);
    }

    #[test]
    fn test_naive_and_fft() {
        let mut coeffs = [FieldElement::new(0);4];
        coeffs[2] = FieldElement::new(5);
        //coeffs[1] = FieldElement::new(1);
        //coeffs[0] = FieldElement::new(11);


        let n = FieldElement::nth_root_of_unity(4);
        let fft =Poly::fft(&coeffs, &n);
        let mut arr = [FieldElement::new(0); 256];
        let len = fft.len();
        arr[..len].copy_from_slice(&fft); // copy the elements from the vector to the array

        let fft2 =Poly::naive_fft(&coeffs, &n);
        let mut arr2 = [FieldElement::new(0); 256];
        let len2 = fft2.len();
        arr2[..len2].copy_from_slice(&fft2); // copy the elements from the vector to the array
        
        assert_eq!(arr[..len], arr2[..len2]);
    }

    #[test]
    fn test_inv_fft() {

        let n = FieldElement::nth_root_of_unity(8);
        let coeffs = [FieldElement::new(1),n, n.pow(2), n.pow(3), n.pow(4),n.pow(5),n.pow(6),n.pow(7)];
        
        let inv_fft =Poly::inv_fft(&coeffs, &n);
        print!("inv_fft: {:?}" , inv_fft);
        let mut arr = [FieldElement::new(0); 256];
        let len = inv_fft.len();
        arr[..len].copy_from_slice(&inv_fft); // copy the elements from the vector to the array
        
        let expected_coeffs = [FieldElement::new(0),FieldElement::new(1), FieldElement::new(0), FieldElement::new(0), FieldElement::new(0), FieldElement::new(0), FieldElement::new(0), FieldElement::new(0)];
        
        assert_eq!(arr[..len], expected_coeffs);
    }
}