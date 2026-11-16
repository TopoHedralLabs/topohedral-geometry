//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

use core::panic;
//{{{ crate imports
//}}}
//{{{ std imports
use std::collections::BTreeSet;
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

pub type Index = u64;

#[derive(Debug)]
pub struct IndexResource {
    next_index: Index,
    free_indices: BTreeSet<Index>,
}

pub const NULL_IDX: Index = Index::MAX;

pub fn is_null(idx: Index) -> bool {
    idx == NULL_IDX
}

impl IndexResource {
    pub fn new() -> Self {
        Self {
            next_index: 0u64,
            free_indices: BTreeSet::new(),
        }
    }

    pub fn peek_index(&self) -> Index {
        if let Some(next_idx) = self.free_indices.first() {
            return *next_idx;
        }
        self.next_index
    }

    pub fn acquire_index(&mut self) -> Index {
        if let Some(next_idx) = self.free_indices.pop_first() {
            return next_idx;
        }
        let next_idx = self.next_index;
        self.next_index += 1;
        next_idx
    }

    pub fn return_index(
        &mut self,
        idx: Index,
    ) {
        assert!(self.in_range(idx));
        self.free_indices.insert(idx);
    }

    pub fn is_free(
        &self,
        idx: Index,
    ) -> bool {
        assert!(self.in_range(idx));
        self.free_indices.contains(&idx)
    }

    pub fn is_used(
        &self,
        idx: Index,
    ) -> bool {
        assert!(self.in_range(idx));
        !self.free_indices.contains(&idx)
    }

    pub fn in_range(
        &self,
        idx: Index,
    ) -> bool {
        idx < self.next_index
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Handle {
    index: Index,
    generation: u32,
}

struct Slot<T> {
    generation: u32,
    value: Option<T>,
}

impl<T> Default for Slot<T> {
    fn default() -> Self {
        Self {
            generation: 0,
            value: None,
        }
    }
}

pub struct SlotVec<T> {
    indices: IndexResource,
    entities: Vec<Slot<T>>,
}

impl<T> SlotVec<T> {
    pub fn new(capacity: usize) -> Self {
        Self {
            indices: IndexResource::new(),
            entities: Vec::<Slot<T>>::with_capacity(capacity),
        }
    }

    pub fn len(&self) -> usize {
        self.entities.len()
    }

    pub fn is_stale(
        &self,
        handle: Handle,
    ) -> bool {
        self.entities[handle.index as usize].generation != handle.generation
    }

    pub fn peek_index(&self) -> Index {
        return self.indices.peek_index();
    }

    pub fn insert(
        &mut self,
        item: T,
    ) -> Handle {
        let idx = self.indices.acquire_index();

        if (idx as usize) < self.entities.len() {
            if self.entities[idx as usize].value.is_none() {
                let generation = self.entities[idx as usize].generation + 1;
                self.entities[idx as usize] = Slot {
                    generation: generation,
                    value: Some(item),
                };
                return Handle {
                    index: idx,
                    generation,
                };
            }
            panic!("Trying to re-take slot");
        }
        if idx as usize == self.entities.len() {
            self.entities.push(Slot {
                generation: 0,
                value: Some(item),
            });
            return Handle {
                index: idx,
                generation: 0,
            };
        }
        panic!("Should not be reachable")
    }

    pub fn get(
        &self,
        handle: Handle,
    ) -> Option<&T> {
        let Handle { index, generation } = handle;
        if self.indices.is_used(index) && self.entities[index as usize].generation == generation {
            return self.entities[index as usize].value.as_ref();
        }
        None
    }

    pub fn get_mut(
        &mut self,
        handle: Handle,
    ) -> Option<&mut T> {
        let Handle { index, generation } = handle;
        if self.indices.is_used(index) && self.entities[index as usize].generation == generation {
            return self.entities[index as usize].value.as_mut();
        }
        None
    }

    pub fn delete(
        &mut self,
        handle: Handle,
    ) {
        let Handle { index, generation } = handle;

        if self.indices.is_free(index) {
            return;
        }
        if generation != self.entities[index as usize].generation {
            return;
        }
        self.entities[index as usize].value = None;
        self.indices.return_index(index);
    }

    pub fn iter(&self) -> Iter<'_, T> {
        let mut start = 0usize;
        while self.indices.is_free(start as Index) {
            start += 1;
        }
        Iter {
            slot_vec: self,
            cur_index: start,
        }
    }
}

pub struct Iter<'a, T> {
    slot_vec: &'a SlotVec<T>,
    cur_index: usize,
}

impl<'a, T> Iterator for Iter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        while self.cur_index < self.slot_vec.len()
            && self.slot_vec.indices.is_free(self.cur_index as Index)
        {
            self.cur_index += 1;
        }

        if self.cur_index == self.slot_vec.len() {
            return None;
        }
        self.slot_vec.entities[self.cur_index].value.as_ref()
    }
}

//-------------------------------------------------------------------------------------------------
//{{{ mod: tests
#[cfg(test)]
mod tests {
    use super::Handle;
    use super::IndexResource;
    use super::SlotVec;

    #[test]
    fn test_indexing() {
        let mut idx_resource = IndexResource::new();

        for i in 0..1000 {
            let next_idx = idx_resource.acquire_index();
            assert_eq!(next_idx, i)
        }

        for i in 0..500 {
            idx_resource.return_index(i);
        }

        for i in 0..500 {
            let next_idx = idx_resource.acquire_index();
            assert_eq!(next_idx, i)
        }

        for i in 0..10 {
            idx_resource.return_index(i);
            assert!(idx_resource.is_free(i))
        }

        assert!(!idx_resource.is_free(1001))
    }

    #[test]
    fn test_slotvec() {
        let mut slotvec = SlotVec::<u64>::new(100);

        for i in 0..200 {
            let handle = slotvec.insert(i);
            assert_eq!(handle.index, i);
            assert_eq!(handle.generation, 0);
            assert!(!slotvec.is_stale(handle))
        }

        for i in 20..50 {
            slotvec.delete(Handle {
                index: i,
                generation: 0,
            });
        }
        assert_eq!(slotvec.peek_index(), 20);

        for i in 20..50 {
            let handle = slotvec.insert(100 + i);
            assert_eq!(handle.index, i);
            assert_eq!(handle.generation, 1);
            assert!(!slotvec.is_stale(handle));
            assert_eq!(slotvec.get(handle), Some(&(100 + i)));
        }
    }
}
//}}}
