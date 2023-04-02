// basic modular arithmetic for field elements

use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, SubAssign, MulAssign, DivAssign};

use std::hash::{Hash, Hasher};


// FieldElement is a wrapper around a u64 that represents a field element
// in the finite field F_p where p = 2^64 - 59
#[derive(Copy, Clone)]
pub struct FieldElement {
    value: u64,
}

impl FieldElement {
    // new creates a new FieldElement from a u64
    pub fn new(value: u64) -> FieldElement {
        FieldElement { value: value % 0xffffffffffffffc5 }
    }

    // add, sub and mul implement the basic arithmetic operations
    pub fn add(&self, other: FieldElement) -> FieldElement {
        FieldElement::new(self.value + other.value)
    }

    pub fn sub(&self, other: FieldElement) -> FieldElement {
        FieldElement::new(self.value + 0xffffffffffffffc5 - other.value)
    }

    pub fn mul(&self, other: FieldElement) -> FieldElement {
        FieldElement::new(self.value * other.value)
    }

    // pow computes the power of a FieldElement
    pub fn pow(&self, exponent: u64) -> FieldElement {
        let mut result = FieldElement::new(1);
        let mut base = self.clone();
        let mut exp = exponent;
        while exp > 0 {
            if exp & 1 == 1 {
                result *= base;
            }
            exp >>= 1;
            base *= base;
        }
        result
    }

    // sqrt computes the square root of a FieldElement
    pub fn sqrt(&self) -> FieldElement {
        self.pow(0x7ffffffffffffffb)
    }
}