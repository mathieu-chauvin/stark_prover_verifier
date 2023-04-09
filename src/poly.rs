// implement polynomials using the FieldElement type

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
            result[i] = self.coeffs[i].add(other.coeffs[i]);
        }
        Poly::new(result)
    }

    fn sub(&self, other: Poly) -> Poly {
        let mut result = [FieldElement::new(0); 256];
        for i in 0..256 {
            result[i] = self.coeffs[i].sub(other.coeffs[i]);
        }
        Poly::new(result)
    }

    //interpolation
    
}