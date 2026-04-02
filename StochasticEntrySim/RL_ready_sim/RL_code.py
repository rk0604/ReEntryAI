import math
import AtmosphereModel
import constants
from interval_math import Interval, promote, box_add, box_scalar_mul, scalar_times_interval
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Protocol, Any, Callable, List
from control import ReentryState

class ReentryEnv:
    def __init__(
        self,
        sim_cfg: SimConfig,
        params: Dict[str, Any],
        control: SigmaControlStack,
    ):
        ...

    def reset(self, initial_state: ReentryState) -> Dict[str, float]:
        ...

    def step(self, action: Any) -> Tuple[Dict[str, float], float, bool, Dict[str, Any]]:
        """
        action can be:
          - sigma_cmd (continuous)
          - delta_sigma_cmd
          - discrete bank reversal
        """
        ...