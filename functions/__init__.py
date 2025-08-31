# Functions package initialization
# Import commonly used classes for easy access

# Process Logic
from .process_logic.cif_processor import CIFProcessor
from .process_logic.supercell_processor import SupercellProcessor
from .process_logic.supercell_verifier import SupercellVerifier
from .process_logic.data_validator import DataValidator
from .process_logic.plane_equation_generator import PlaneEquationGenerator
from .process_logic.atom_plane_assignment import AtomPlaneAssignment
from .process_logic.reference_plane_selector import ReferencePlaneSelector
from .process_logic.plane_comparison_analyzer import PlaneComparisonAnalyzer
from .process_logic.reference_plane_manager import ReferencePlaneManager
from .process_logic.supercell_expander import SupercellExpander
from .process_logic.miller_verification_manager import MillerVerificationManager
from .process_logic.sequence_and_cell_calculator import SequenceAndCellCalculator
from .process_logic.qe_output_generator import QEOutputGenerator
from .process_logic.vacuum_processor import VacuumProcessor
from .process_logic.coordinate_transformer import CoordinateTransformer
from .process_logic.structure_processor import StructureProcessor
from .process_logic.hydrogen_placer import HydrogenPlacer
from .process_logic.hydrogen_optimizer import HydrogenOptimizer

# Process Logic - Writers
from .process_logic.writers.qe_writer import QEWriter
from .process_logic.writers.supercell_writer import SupercellWriter
from .process_logic.writers.plane_writer import PlaneWriter
from .output_logic.display_manager_base import BaseDisplayManager
from .output_logic.display_manager_phase1 import Phase1DisplayManager
from .output_logic.display_manager_phase2 import Phase2DisplayManager
from .output_logic.display_manager_phase3 import Phase3DisplayManager
from .output_logic.display_manager_phase4 import Phase4DisplayManager
from .output_logic.display_manager_phase5 import Phase5DisplayManager
from .output_logic.console_writer import ConsoleWriter
from .utils.temporary_file_handler import TemporaryFileHandler

# Input Logic
from .input_logic.cif_reader import CIFReader
from .input_logic.data_parser import CIFParser
from .input_logic.miller_index_input import MillerIndexInput
from .input_logic.vacuum_layer_input import VacuumLayerInput
from .input_logic.hydrogenation_input import HydrogenationInput

# Utils
from .utils.coordinates import convert_atoms_to_cartesian

# Factory
from .factory import ComponentFactory
