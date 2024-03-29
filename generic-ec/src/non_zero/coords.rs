use crate::coords::{
    AlwaysHasAffineX, AlwaysHasAffineXY, AlwaysHasAffineY, HasAffineX, HasAffineXY, HasAffineY,
};
use crate::{Curve, Point};

use super::definition::NonZero;

impl<E: Curve> AlwaysHasAffineX<E> for NonZero<Point<E>>
where
    Point<E>: HasAffineX<E>,
{
    fn x(&self) -> crate::coords::Coordinate<E> {
        // The only point that may not have coords for some curves is point at infinity (or
        // identity point). Since we know it's non-zero, it must have coordinates
        #![allow(clippy::expect_used)]
        HasAffineX::x(&**self).expect("non-zero point always has coordinates")
    }
}

impl<E: Curve> AlwaysHasAffineY<E> for NonZero<Point<E>>
where
    Point<E>: HasAffineY<E>,
{
    fn y(&self) -> crate::coords::Coordinate<E> {
        // The only point that may not have coords for some curves is point at infinity (or
        // identity point). Since we know it's non-zero, it must have coordinates
        #![allow(clippy::expect_used)]
        HasAffineY::y(&**self).expect("non-zero point always has coordinates")
    }
}

impl<E: Curve> AlwaysHasAffineXY<E> for NonZero<Point<E>>
where
    Point<E>: HasAffineXY<E>,
{
    fn from_coords(coords: &crate::coords::Coordinates<E>) -> Option<Self> {
        <Point<E> as HasAffineXY<E>>::from_coords(coords).and_then(NonZero::from_point)
    }
}
