// basic modular arithmetic for field elements

use std::env::temp_dir;
use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, SubAssign, MulAssign, DivAssign};

use std::hash::{Hash, Hasher};


    // FieldElement is a wrapper around a u64 that represents a field element
    // in the finite field F_p where p = 2^64 - 2^32 +1


//p in hex
// 2^64 - 2^32 +1
// 2^64 = 0x10000000000000000
// 2^32 = 0x100000000
// 2^64 - 2^32 = 0x10000000000000000 - 0x100000000 = 0xffffffff00000000
// 2^64 - 2^32 +1 = 0xffffffff00000000 + 1 = 0xffffffff00000001

// generator = 2


// declare a constant for p
const P: u64 = 0xffffffff00000001;

#[derive(Copy, Clone, Debug)]
pub struct FieldElement {
    value: u64,
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &FieldElement) -> bool {
        self.value == other.value
    }
}

impl Add<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, other: FieldElement) -> FieldElement {
        FieldElement::new(self.value + other.value)
    }
}

impl Sub<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn sub(self, other: FieldElement) -> FieldElement {
        FieldElement::new(self.value + P - other.value)
    }
}

impl Mul<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn mul(self, other: FieldElement) -> FieldElement {
        let res = (self.value as u128 * other.value as u128) % P as u128 ;
        FieldElement::new(res.try_into().unwrap())
    }
}

impl FieldElement {
    // new creates a new FieldElement from a u64
    pub fn new(value: u64) -> FieldElement {
        FieldElement { value: value % P }
    }

    // pow computes the exponentiation of a FieldElement by using the binary exponentiation algorithm
    pub fn pow(&self, exp: u64) -> FieldElement {
        let mut result = FieldElement::new(1);
        let mut base = self.value;
        let mut exp = exp;
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * FieldElement::new(base);
            }
            exp = exp >> 1;
            base = (base * base) % P;
        }
        result
    }
    

    // inv computes the multiplicative inverse of a FieldElement by using the extended Euclidean algorithm
    pub fn inv(&self) -> FieldElement {
        let mut r0 = self.value% P;
        let mut r1 = P;
        let mut s0:i128 = 1;
        let mut s1:i128 = 0;
        
        let mut q;
        let mut tmp;
        let mut tmp2:i128;
        while r0 > 1 {
            q = r1 / r0;

            //println!("r0: {}, r1: {}, s0: {}, s1: {}, q: {}", r0, r1, s0, s1, q);

            tmp = r1;
            r1 = r0;
            r0 = tmp - q * r0;

            tmp2 = s1;
            s1 = s0;
            s0 = tmp2 - q as i128 * s0;
            
        }

        //find s0 modulo P
        if s0 < 0 {
            s0 = s0 + P as i128;
        }

        return FieldElement::new(s0.try_into().unwrap());

    }

    // find multiples inverses using montgomery batch inversion
    fn multi_inv(&self, values: &[FieldElement]) -> Vec<FieldElement> {
        // declare a empty vector to hold the partials
        let mut partials = Vec::new();
        // set the first partial to a
        partials.push(values[0]);
        // for i in range (0,values) 
        for i in 1..values.len() {
            // get the value at index i
            let value = values[i];
            // compute the partial product
            let partial:FieldElement = partials[i-1] * value;
            // append the partial product to the partials vector
            partials.push(partial);
        }


        // calculate the inverse of the last partial product
        let mut inv = partials[partials.len()-1].inv();

        // calculate outputs

        // initiate an empty vector with an allocate size of partial len values

        let mut outputs = Vec::new();

        
        for i in (1..values.len()).rev() {
            print!("i: {}", i);
            let output = inv * partials[i-1];
            print!("output: {}", output.value);
            // set the ith element of outputs to output
            outputs.push(output);
            inv = inv * values[i];
        }

        // set the first element of outputs to inv
        outputs.push(inv);

        // reverse the outputs vector
        outputs.reverse();

        outputs

        
    }
    

    
}



#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_add() {
        let a = FieldElement::new(1);
        let b = FieldElement::new(2);
        assert_eq!(a+b, FieldElement::new(3));
    }

    #[test]
    fn test_add_overflow() {
        let a = FieldElement::new(P-1);
        let b = FieldElement::new(2);
        assert_eq!(a+b, FieldElement::new(1));
    }


    #[test]
    fn test_sub() {
        let a = FieldElement::new(3);
        let b = FieldElement::new(2);
        assert_eq!(a-b, FieldElement::new(1));
    }

    #[test]
    fn test_sub_overflow() {
        let a = FieldElement::new(1);
        let b = FieldElement::new(2);
        assert_eq!(a-b, FieldElement::new(P-1));
    }

    #[test]
    fn test_inv() {
        let a = FieldElement::new(42);
        let inv = a.inv();
        assert_eq!(a*inv, FieldElement::new(1));
    }

    #[test]
    fn test_pow() {
        let a = FieldElement::new(2);
        let b = 3 as u64;
        assert_eq!(a.pow(b), FieldElement::new(8));
    }

    #[test]
    fn test_multi_inv() {
        let a = FieldElement::new(42);
        let b = FieldElement::new(17);
        let c = FieldElement::new(13);
        let invs = a.multi_inv(&[a,b,c]);
        assert_eq!(a*invs[0], FieldElement::new(1));
        assert_eq!(b*invs[1], FieldElement::new(1));
        assert_eq!(c*invs[2], FieldElement::new(1));
    }

}



