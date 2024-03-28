use core::cmp;
use core::hash::{self, Hash};

use elliptic_curve::group::cofactor::CofactorGroup;
use elliptic_curve::{
    sec1::{EncodedPoint, FromEncodedPoint, ModulusSize, ToEncodedPoint},
    CurveArithmetic, FieldBytesSize, Group,
};
use generic_ec_core::*;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use zeroize::Zeroize;

pub struct RustCryptoPoint<E: CurveArithmetic>(pub E::ProjectivePoint);

impl<E: CurveArithmetic> Additive for RustCryptoPoint<E> {
    #[inline]
    fn add(a: &Self, b: &Self) -> Self {
        RustCryptoPoint(a.0 + b.0)
    }

    #[inline]
    fn sub(a: &Self, b: &Self) -> Self {
        RustCryptoPoint(a.0 - b.0)
    }

    #[inline]
    fn negate(x: &Self) -> Self {
        RustCryptoPoint(-x.0)
    }

    #[inline]
    fn double(x: &Self) -> Self {
        RustCryptoPoint(x.0.double())
    }
}

impl<E: CurveArithmetic> From<CurveGenerator> for RustCryptoPoint<E> {
    #[inline]
    fn from(_: CurveGenerator) -> Self {
        RustCryptoPoint(E::ProjectivePoint::generator())
    }
}

impl<E: CurveArithmetic> Zero for RustCryptoPoint<E> {
    #[inline]
    fn zero() -> Self {
        RustCryptoPoint(E::ProjectivePoint::identity())
    }

    #[inline]
    fn is_zero(x: &Self) -> Choice {
        x.0.is_identity()
    }
}

impl<E: CurveArithmetic> OnCurve for RustCryptoPoint<E> {
    #[inline]
    fn is_on_curve(&self) -> Choice {
        Choice::from(1)
    }
}

impl<E> SmallFactor for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: CofactorGroup,
{
    #[inline]
    fn is_torsion_free(&self) -> Choice {
        self.0.is_torsion_free()
    }
}

impl<E> ConstantTimeEq for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: ConstantTimeEq,
{
    #[inline]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl<E> ConditionallySelectable for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: ConditionallySelectable,
{
    #[inline]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self(E::ProjectivePoint::conditional_select(&a.0, &b.0, choice))
    }

    #[inline]
    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        self.0.conditional_assign(&other.0, choice)
    }

    #[inline]
    fn conditional_swap(a: &mut Self, b: &mut Self, choice: Choice) {
        E::ProjectivePoint::conditional_swap(&mut a.0, &mut b.0, choice)
    }
}

impl<E> CompressedEncoding for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::AffinePoint: ToEncodedPoint<E> + From<E::ProjectivePoint>,
    FieldBytesSize<E>: ModulusSize,
{
    type Bytes = elliptic_curve::sec1::CompressedPoint<E>;
    fn to_bytes_compressed(&self) -> Self::Bytes {
        let point_encoded = E::AffinePoint::from(self.0).to_encoded_point(true);

        let mut bytes = Self::Bytes::default();
        if !bool::from(Self::is_zero(self)) {
            bytes.copy_from_slice(point_encoded.as_bytes());
        }

        bytes
    }
}

impl<E> UncompressedEncoding for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::AffinePoint: ToEncodedPoint<E> + From<E::ProjectivePoint>,
    FieldBytesSize<E>: ModulusSize,
{
    type Bytes = elliptic_curve::sec1::UncompressedPoint<E>;
    fn to_bytes_uncompressed(&self) -> Self::Bytes {
        let point_encoded = E::AffinePoint::from(self.0).to_encoded_point(false);

        let mut bytes = Self::Bytes::default();
        if !bool::from(Self::is_zero(self)) {
            bytes.copy_from_slice(point_encoded.as_bytes());
        }

        bytes
    }
}

impl<E> Decode for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::AffinePoint: FromEncodedPoint<E> + Into<E::ProjectivePoint>,
    FieldBytesSize<E>: ModulusSize,
{
    fn decode(mut bytes: &[u8]) -> Option<Self> {
        let all_zero = bytes.iter().all(|b| *b == 0);
        if all_zero {
            // This is the only representation of identity point recognized
            // by `elliptic-curve` library
            bytes = &[0]
        }
        let encoded_point = EncodedPoint::<E>::from_bytes(bytes).ok()?;
        Option::from(E::AffinePoint::from_encoded_point(&encoded_point))
            .map(|point: E::AffinePoint| Self(point.into()))
    }
}

impl<E> Clone for RustCryptoPoint<E>
where
    E: CurveArithmetic,
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<E> Copy for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: Copy,
{
}

impl<E> Zeroize for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: Zeroize,
{
    fn zeroize(&mut self) {
        self.0.zeroize()
    }
}

impl<E> PartialEq for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<E> Eq for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: Eq,
{
}

impl<E> Hash for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: Hash,
{
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state)
    }
}

impl<E> PartialOrd for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<E> Ord for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: Ord,
{
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl<E> Default for RustCryptoPoint<E>
where
    E: CurveArithmetic,
    E::ProjectivePoint: Default,
{
    fn default() -> Self {
        Self(Default::default())
    }
}
