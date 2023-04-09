// basic modular arithmetic for field elements

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


impl FieldElement {
    // new creates a new FieldElement from a u64
    pub fn new(value: u64) -> FieldElement {
        FieldElement { value: value % P }
    }


    

    // inv computes the multiplicative inverse of a FieldElement by using the extended Euclidean algorithm
    pub fn inv(&self) -> FieldElement {
        let mut u1 = 1u64;
        let mut u2 = 0u64;
        let mut u3 = self.value;
        let mut v1 = 0u64;
        let mut v2 = 1u64;
        let mut v3 = Pu64;
        while v3 != 0 {
            let q = u3 / v3;
            let t1 = u1 - q * v1;
            let t2 = u2 - q * v2;
            let t3 = u3 - q * v3;
            u1 = v1;
            u2 = v2;
            u3 = v3;
            v1 = t1;
            v2 = t2;
            v3 = t3;
        }
        FieldElement::new(u1)
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

}



