from pathlib import Path
from typing import Optional, Union
from pydantic import BaseModel, Field

PathLike = Union[str, Path]


class FluidPhysicalProperties(BaseModel):
    rho: float
    nu: float

FLUID_PYSICAL_PROPERTIES = {
    "Air": FluidPhysicalProperties(rho= 1.225, nu= 1.5e-5),
    "Water": FluidPhysicalProperties(rho= 1000, nu= 1e-6),
}

class StlInputModel(BaseModel):
    stl_path: PathLike
    purpose: str

class ProjectInputModel(BaseModel):
    project_path: PathLike
    refinement_level: int
    is_internal_flow: bool
    on_ground: Optional[bool] = None
    fluid_properties: FluidPhysicalProperties
    inlet_values: tuple[float, float, float]
    is_transient: bool
    n_core: int
    is_half_model: bool
    use_function_objects: bool
    stl_files: list[StlInputModel] = Field(default_factory=list)


