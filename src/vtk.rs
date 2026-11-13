//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::Vec3;
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use vtkio::model::{
    Attribute, ByteOrder, CellType, Cells, DataSet, UnstructuredGridPiece, Version, VertexNumbers,
    Vtk,
};
//}}}
//--------------------------------------------------------------------------------------------------

pub(crate) trait ToVtk {
    type Options;
    fn add_to_vtk(
        &self,
        opts: Self::Options,
        vtk_file: &mut VtkUgridFile,
    );
}

#[derive(Default)]
pub(crate) struct VtkUgridFile {
    coords: Vec<f64>,
    connectivity: Vec<u64>,
    offsets: Vec<u64>,
    types: Vec<CellType>,
    point_data: Vec<Attribute>,
    cell_data: Vec<Attribute>,
}

impl VtkUgridFile {
    /// Running point count: callers use this as the base index for their cells, exactly as
    /// `Surface::save` does with `np_init{vtk_file.numPoints()}`.
    pub(crate) fn num_points(&self) -> usize {
        self.coords.len() / 3
    }

    pub(crate) fn num_cells(&self) -> usize {
        self.types.len()
    }

    pub(crate) fn add_point(
        &mut self,
        p: Vec3,
    ) {
        self.coords.extend_from_slice(p.as_ref());
    }

    pub(crate) fn add_cell(
        &mut self,
        cell: Vec<u64>,
        cell_type: CellType,
    ) {
        assert!(
            cell.iter().all(|&i| (i as usize) < self.num_points()),
            "point id out of range"
        );
        self.connectivity.extend(cell);
        self.offsets.push(self.connectivity.len() as u64);
        self.types.push(cell_type);
    }

    pub(crate) fn save(
        self,
        path: &str,
    ) -> Result<(), vtkio::Error> {
        Vtk {
            version: Version::new((0, 1)),
            byte_order: ByteOrder::LittleEndian,
            title: String::from("topohedral"),
            file_path: None,
            data: DataSet::inline(UnstructuredGridPiece {
                points: self.coords.into(),
                cells: Cells {
                    cell_verts: VertexNumbers::XML {
                        connectivity: self.connectivity,
                        offsets: self.offsets,
                    },
                    types: self.types,
                },
                data: vtkio::model::Attributes {
                    point: self.point_data,
                    cell: self.cell_data,
                },
            }),
        }
        .export(path)
    }
}
