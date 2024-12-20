from pathlib import Path
import shutil
from typing import Union
from src.primitives import AmpersandIO
from src.project import AmpersandProject


class ProjectService:



    @staticmethod
    def create_project(project_path: Union[str, Path]):
        """Create the project directory structure for an OpenFOAM case.
        
        Args:
            project_path: Path where the project should be created, can be string or Path

        Raises:
            ValueError: If project_path is None
            OSError: If there are issues creating directories or changing directory
        """
        project = AmpersandProject(project_path)

        AmpersandIO.printMessage("Creating the project")

        if project_path is None:
            raise ValueError("No project path provided")

        project_path = Path(project_path)

        # Create required OpenFOAM directories 
        required_dirs = [
            project_path / "0",
            project_path / "constant",
            project_path / "system",
            project_path / "constant/triSurface"
        ]

        try:
            for directory in required_dirs:
                if not directory.exists():
                    directory.mkdir(parents=True)
                    AmpersandIO.printMessage(f"Created {directory} directory")
        except OSError as e:
            raise OSError(f"Failed to create OpenFOAM directory structure: {e}")
        
        project.load_default_settings()
        project.write_settings()

        return project

    @staticmethod
    def load_project(project_path: Union[str, Path]):
        project_path = Path(project_path)
        project = AmpersandProject(project_path)

        AmpersandIO.printMessage(f"Project path: {project_path}")
        if not project_path.exists():
            raise FileNotFoundError("No project found. Exiting the program")
        AmpersandIO.printMessage("Loading the project")

        project.load_settings()
        ProjectService.validate_project(project)

        AmpersandIO.printMessage("Project loaded successfully")

        return project


    @staticmethod
    def validate_project(project: Union[AmpersandProject, str, Path]):
        project_path = Path(project) if isinstance(project, (str, Path)) else project.project_path
        ProjectService.check_directory(project_path / "0", project_path / "0.orig", copy=True)
        ProjectService.check_directory(project_path / "constant")
        ProjectService.check_directory(project_path / "system")
        ProjectService.check_directory(project_path / "constant/triSurface", check_files=True)



    @staticmethod
    def check_directory(path: Union[Path, str], fallback=None, copy=False, check_files=False):
        path = Path(path)
        if not path.exists():
            if fallback and Path(fallback).exists() and copy:
                AmpersandIO.printMessage(f"{fallback} directory found. Copying to {path}")
                shutil.copytree(fallback, path)
            else:
                raise FileNotFoundError(f"{path} directory not found.")
        if check_files:
            if not list(path.iterdir()):
                raise FileNotFoundError(f"No files found in {path}.")

    @staticmethod
    def check_log_files(project_path: Union[str, Path]):
        project_path = Path(project_path)
        log_files = list(project_path.iterdir())
        if 'log.simpleFoam' in [f.name for f in log_files] or 'log.pimpleFoam' in [f.name for f in log_files]:
            AmpersandIO.printMessage("Simulation log file found")
            return True
        AmpersandIO.printMessage("No simulation log files found.")
        return False

    @staticmethod
    def check_post_process_files(project_path: Union[str, Path]):
        project_path = Path(project_path)
        post_process_path = project_path / "postProcessing/probe/0"
        if not post_process_path.exists():
            AmpersandIO.printMessage(f"{post_process_path} directory does not exist.")
            return False
        files = [f.name for f in post_process_path.iterdir()]
        if 'U' not in files or 'p' not in files:
            AmpersandIO.printMessage("Required files 'U' and 'p' not found in postProcessing.")
            return False
        return True

    @staticmethod
    def check_forces_files(project_path: Union[str, Path]):
        project_path = Path(project_path)
        forces_path = project_path / "postProcessing/forces/0"
        if not forces_path.exists():
            AmpersandIO.printMessage(f"{forces_path} directory does not exist.")
            return False
        files = [f.name for f in forces_path.iterdir()]
        if 'force.dat' not in files:
            AmpersandIO.printMessage("force.dat file not found in forces directory.")
            return False
        return True
