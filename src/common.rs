//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
//}}}
//{{{ std imports
use std::{
    num,
    ops::{Index, IndexMut},
};
//}}}
//{{{ dep imports
use topohedral_linalg::SCVector;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ types: Vec3 and Vec4
/// A point or vector in Euclidean 3-space.
pub type Vec3 = SCVector<f64, 3>;
/// A point or vector in homogeneous 3-space.
pub type Vec4 = SCVector<f64, 4>;
//}}}
//{{{ fun: inv_homog
/// Performs the perspective map (inverse of homogeneuos map) from homogeneous coordinates to
/// Euclidean coordinates.
pub(crate) fn inv_homog(point_w: &Vec4) -> Vec3 {
    let mut point = Vec3::zeros();
    let w = point_w[3];
    for i in 0..3 {
        point[i] = point_w[i] / w;
    }
    point
}
//}}}
//{{{ fun: homog
/// Performs the inverse perspective map (homogeneous map) from Euclidean coordinates to
/// Homogeneious coordinates.
pub(crate) fn homog(
    point: &Vec3,
    weight: f64,
) -> Vec4 {
    let mut point_w = Vec4::zeros();
    for i in 0..3 {
        point_w[i] = weight * point[i];
    }
    point_w[3] = weight;
    point_w
}
//}}}
//{{{ struct: IndexHelper
/// This struct provides the index-conversion for going between a 1D index and an N-dimensional index.
///
/// Indexing in this struct is implemented such that the leftmost index varies the fastest.
/// In 2D this is equivalent to column-major ordering. So, given a set of dimensions:
/// $$
/// (n_{1}, n_{2}, ..., n_{N})
/// $$
/// The linear index $j$ of a multi-dimensional index $(i_{1}, i_{2}, ..., i_{N})$ is given by:
/// $$
/// j = i_{1} + n_{1} i_{2} + n_{1}n_{2}i_{3} + ... + n_{1}n_{2}...n_{N-1}i_{N}
/// $$
/// And given a linear index $j$, the multi-dimensional index $(i_{1}, i_{2}, ..., i_{N})$ is:
/// $$
/// todo
/// $$
#[repr(C)]
pub(crate) struct IndexHelper<const N: usize> {
    /// Helper array for converting a multi-dimensional index into a linear index
    lin_helper: [usize; N],
}
//}}}
//{{{ impl: IndexHelper
impl<const N: usize> IndexHelper<N> {
    pub fn new(dims: &[usize]) -> Self {
        debug_assert!(dims.len() >= N);

        let mut helper = IndexHelper::<N> { lin_helper: [1; N] };

        for i in 1..N {
            helper.lin_helper[i] = helper.lin_helper[i - 1] * dims[i - 1];
        }
        helper
    }

    pub fn lin_index(
        &self,
        indices: &[usize],
    ) -> usize {
        debug_assert!(indices.len() >= N);
        let mut idx = 0usize;
        for (&index, &stride) in indices.iter().zip(&self.lin_helper).take(N) {
            idx += index * stride;
        }
        idx
    }
}
//}}}
//{{{ struct: StaticArray
pub(crate) struct StatArray<T, const MAX_SIZE: usize, const DIM: usize>
where
    T: Clone,
{
    data: [T; MAX_SIZE],
    index_helper: IndexHelper<DIM>,
}
//}}}
//{{{ impl: StaticArray
impl<T, const MAX_SIZE: usize, const DIM: usize> StatArray<T, MAX_SIZE, DIM>
where
    T: Clone,
{
    pub fn new(
        init_value: T,
        dims: &[usize],
    ) -> Self {
        let num_elems: usize = dims.iter().product();
        assert!(dims.len() == DIM);
        assert!(num_elems <= MAX_SIZE);

        Self {
            data: std::array::from_fn(|_| init_value.clone()),
            index_helper: IndexHelper::new(dims),
        }
    }

    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    pub fn as_slice_mut(&mut self) -> &mut [T] {
        &mut self.data
    }
}
//}}}
//{{{ impl: Index for StaticArray
impl<T, const MAX_SIZE: usize, const DIM: usize> Index<&[usize; DIM]>
    for StatArray<T, MAX_SIZE, DIM>
where
    T: Clone,
{
    type Output = T;

    fn index(
        &self,
        index: &[usize; DIM],
    ) -> &Self::Output {
        let lin_index = self.index_helper.lin_index(index);
        &self.data[lin_index]
    }
}
//}}}
//{{{ impl: IndexMut for StaticArray
impl<T, const MAX_SIZE: usize, const DIM: usize> IndexMut<&[usize; DIM]>
    for StatArray<T, MAX_SIZE, DIM>
where
    T: Clone,
{
    fn index_mut(
        &mut self,
        index: &[usize; DIM],
    ) -> &mut Self::Output {
        let lin_index = self.index_helper.lin_index(index);
        &mut self.data[lin_index]
    }
}
//}}}
//{{{ struct: Array
pub(crate) struct DynArray<T, const DIM: usize>
where
    T: Clone,
{
    data: Vec<T>,
    index_helper: IndexHelper<DIM>,
}
//}}}
//{{{ impl: Array
impl<T, const DIM: usize> DynArray<T, DIM>
where
    T: Clone,
{
    pub fn new(
        init_value: T,
        dims: &[usize],
    ) -> Self {
        let num_elems: usize = dims.iter().product();
        assert!(dims.len() == DIM);

        Self {
            data: vec![init_value; num_elems],
            index_helper: IndexHelper::new(dims),
        }
    }

    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    pub fn as_slice_mut(&mut self) -> &mut [T] {
        &mut self.data
    }
}
//}}}
//{{{ impl: Index for Array
impl<T, const DIM: usize> Index<&[usize; DIM]> for DynArray<T, DIM>
where
    T: Clone,
{
    type Output = T;

    fn index(
        &self,
        index: &[usize; DIM],
    ) -> &Self::Output {
        let lin_index = self.index_helper.lin_index(index);
        &self.data[lin_index]
    }
}
//}}}
//{{{ impl: IndexMut for Array
impl<T, const DIM: usize> IndexMut<&[usize; DIM]> for DynArray<T, DIM>
where
    T: Clone,
{
    fn index_mut(
        &mut self,
        index: &[usize; DIM],
    ) -> &mut Self::Output {
        let lin_index = self.index_helper.lin_index(index);
        &mut self.data[lin_index]
    }
}
//}}}
//{{{ fun: binom_coeff
pub(crate) fn binom_coeff<const MAX_SIZE: usize>(n: usize) -> StatArray<f64, MAX_SIZE, 2> {
    let mut binom = StatArray::<f64, MAX_SIZE, 2>::new(0.0, &[n + 1, n + 1]);
    for i in 0..n + 1 {
        binom[&[i, i]] = 1.0;
        binom[&[i, 0]] = 1.0;
    }
    for n2 in 2..n + 1 {
        for k2 in 1..n2 {
            binom[&[n2, k2]] = binom[&[n2 - 1, k2 - 1]] + binom[&[n2 - 1, k2]];
        }
    }
    binom
}
//}}}
