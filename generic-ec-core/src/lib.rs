//! # Core traits of `generic-ec` crate
//!
//! This crate is not currently properly documented, and API is not considered stable.

#![no_std]
#![cfg_attr(not(test), forbid(unused_crate_dependencies))]
#![cfg_attr(not(test), deny(clippy::unwrap_used, clippy::expect_used))]

use core::fmt::Debug;
use core::hash::Hash;

use generic_array::{ArrayLength, GenericArray};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};
use zeroize::Zeroize;

pub mod coords;
pub mod hash_to_curve;

/// Elliptic curve
///
/// This trait contains all the low-level curve implementation logic: scalar, point arithmetics,
/// encoding and etc.
pub trait Curve: Debug + Copy + Eq + Ord + Hash + Default + Sync + Send + 'static {
    const CURVE_NAME: &'static str;

    type Point: Additive
        + From<CurveGenerator>
        + Zero
        + Zeroize
        + OnCurve
        + SmallFactor
        + Copy
        + Eq
        + ConstantTimeEq
        + ConditionallySelectable
        + Default
        + CompressedEncoding<Bytes = Self::CompressedPointArray>
        + UncompressedEncoding<Bytes = Self::UncompressedPointArray>
        + Decode
        + Unpin
        + Sync
        + Send;
    type Scalar: Additive
        + Multiplicative<Self::Scalar, Output = Self::Scalar>
        + Multiplicative<CurveGenerator, Output = Self::Point>
        + Multiplicative<Self::Point, Output = Self::Point>
        + Invertible
        + Zero
        + One
        + Samplable
        + Zeroize
        + Copy
        + Eq
        + ConstantTimeEq
        + ConditionallySelectable
        + Default
        + IntegerEncoding<Bytes = Self::ScalarArray>
        + Unpin
        + Sync
        + Send;

    type CompressedPointArray: ByteArray;
    type UncompressedPointArray: ByteArray;
    type ScalarArray: ByteArray;
    type CoordinateArray: ByteArray;
}

pub trait Additive {
    fn add(a: &Self, b: &Self) -> Self;
    fn sub(a: &Self, b: &Self) -> Self;
    fn negate(x: &Self) -> Self;

    /// Takes `x`, returns `x + x`
    ///
    /// This can be more efficient than calling [`Self::add(x, x)`](Self::add)
    fn double(x: &Self) -> Self
    where
        Self: Sized,
    {
        Self::add(x, x)
    }
}

pub trait Multiplicative<Rhs> {
    type Output;
    fn mul(a: &Self, b: &Rhs) -> Self::Output;
}

pub trait Invertible
where
    Self: Sized,
{
    fn invert(x: &Self) -> CtOption<Self>;
}

/// Type that has zero value (additive identity)
pub trait Zero {
    /// Constructs zero value of `Self`
    fn zero() -> Self;
    /// Checks (in constant-time) if `x` is zero
    fn is_zero(x: &Self) -> Choice;
}

/// Type that has "one" value (multiplicative identity)
pub trait One {
    /// Constructs one value of `Self`
    fn one() -> Self;
    /// Checks (in constant-time) if `x` is one
    fn is_one(x: &Self) -> Choice;
}

/// Type can be uniformely sampled from source of randomness
pub trait Samplable {
    /// Uniformely samples a random value of `Self`
    fn random<R: RngCore>(rng: &mut R) -> Self;
}

pub trait OnCurve {
    fn is_on_curve(&self) -> Choice;
}

pub trait SmallFactor {
    fn is_torsion_free(&self) -> Choice;
}

pub struct CurveGenerator;

pub trait CompressedEncoding
where
    Self: Sized,
{
    type Bytes: ByteArray;

    fn to_bytes_compressed(&self) -> Self::Bytes;
}

pub trait UncompressedEncoding
where
    Self: Sized,
{
    type Bytes: ByteArray;

    fn to_bytes_uncompressed(&self) -> Self::Bytes;
}

pub trait IntegerEncoding
where
    Self: Sized,
{
    type Bytes: ByteArray;

    fn to_be_bytes(&self) -> Self::Bytes;
    fn to_le_bytes(&self) -> Self::Bytes;

    fn from_be_bytes(bytes: &Self::Bytes) -> Self;
    fn from_le_bytes(bytes: &Self::Bytes) -> Self;

    fn from_be_bytes_exact(bytes: &Self::Bytes) -> Option<Self>;
    fn from_le_bytes_exact(bytes: &Self::Bytes) -> Option<Self>;
}

pub trait Decode: Sized {
    fn decode(bytes: &[u8]) -> Option<Self>;
}

pub struct Error;

pub trait ByteArray: AsRef<[u8]> + AsMut<[u8]> + Clone + Send + Sync + 'static {
    /// New byte array of zeroes
    ///
    /// Alternative to [`Default`] that is not implemented for generic `[T; N]`
    /// (see [tracking issue](https://github.com/rust-lang/rust/issues/61415))
    fn zeroes() -> Self;
}

impl<const N: usize> ByteArray for [u8; N] {
    fn zeroes() -> Self {
        [0; N]
    }
}

impl<N: ArrayLength<u8>> ByteArray for GenericArray<u8, N> {
    fn zeroes() -> Self {
        GenericArray::default()
    }
}
