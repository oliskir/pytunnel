# ================================================================================= #
#    License: GNU GPLv3                                                             #
#                                                                                   #
#        This program is free software: you can redistribute it and/or modify       # 
#        it under the terms of the GNU General Public License as published by       #
#        the Free Software Foundation, either version 3 of the License, or          #
#        (at your option) any later version.                                        #
#                                                                                   #
#        This program is distributed in the hope that it will be useful,            #
#        but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#        GNU General Public License for more details.                               #
#                                                                                   #
#        You should have received a copy of the GNU General Public License          #
#        along with this program.  If not, see <https://www.gnu.org/licenses/>.     #
#                                                                                   #
# ================================================================================= #

""" barriers module

    This module provides implementations of various potential barriers.

    Barriers are implemented as classes, following a standard template.

    In particular, the following methods must be present
    
        has_changed
        eval
        min
        max
        left_edge
        right_edge
        
    and the constructor may only accept a single dictionary as 
    argument. For an example, see the Square barrier class.
        
    Contents:
        Square class: 
        BinnedBarrier class

    Author: Oliver Kirsebom
    Contact: oliver.kirsebom@gmail.com
"""

import math
import numpy as np
import pytunnel.units as units
from pint import UnitRegistry


class SquareStep():
    """ Square barrier with option to suddenly change the barrier height.
        The potential to the left of the barrier can be different from zero.

        Args:
            params: dict
                Barrier parameters:
                    * width: width in fm
                    * initial_height: initial height in MeV  
                    * final_height: final height in MeV 
                    * transition_pos: wave packet classical position in fm when barrier changes height
                    * left_pot: potential left of barrier in MeV
                The first two are mandatory, the last two are optional.
            x_start: float 
                Wave packet start position in fm
            velocity: float
                Wave packet velocity in fm/s 

        Attributes:
            initial_height: float
                Initial height of the barrier in reduced units
            width: float
                Width of the barrier in reduced units
            transition_time: float
                Time at which barrier height changes
            final_height: float
                Final height of the barrier in reduced units
            left_pot: float
                potential left of barrier in reduced units
    """
    def __init__(self, params, x_start, velocity):

        ureg = UnitRegistry()
        Q = ureg.Quantity

        for s in ['initial_height', 'height']:
            if s in params.keys():
                self.initial_height = Q(params[s]).m_as("MeV") / units.MeV
        
        self.width = Q(params['width']).m_as("fm") / units.fm

        if 'transition_pos' in params.keys():
            self.transition_time = (Q(params['transition_pos']).m_as("fm") - x_start) / velocity / units.sec
        else:
            self.transition_time = math.inf

        if 'final_height' in params.keys():
            self.final_height = Q(params['final_height']).m_as("MeV") / units.MeV
        else:
            self.final_height = self.initial_height

        if 'left_pot' in params.keys():
            self.left_pot = Q(params['left_pot']).m_as("MeV") / units.MeV
        else:
            self.left_pot = 0

    def has_changed(self, t1, t2):
        """ Determine if the barrier shape at time t2 has changed 
            compared to the shape at time t1.

            Times are given in reduced units. Multiply by units.sec 
            to convert to seconds.

            Args:
                t1: float
                    Reference time
                t2: float
                    Current time

            Returns:
                res: bool
                    True, if the barrier shape has changed. False, otherwise
        """
        if t1 is None or t2 is None:
            return False

        h1 = self._height(t1)
        h2 = self._height(t2)
        res = (h1 != h2)

        return res

    def eval(self, x, t):
        """ Evaluate the potential barrier at the specified position and time.

            Position and time are given in reduced units. Multiply by 
            units.fm and units.sec to convert to fm and seconds.

            Args:
                x: float or array-like
                    Position
                t: float
                    Time

            Returns:
                val: float or array-like
                    Value(s) of the potential
        """
        h = self._height(t)

        val = np.logical_and(x >= 0, x < self.width)
        val = val.astype(np.float)
        val *= h

        if self.left_pot != 0: val[x < 0] = self.left_pot

        return val

    def min(self):
        """ Minimum value of the potential barrier for all positions and times.

            Returns:
                v: float
                    Minimum barrier height in reduced units
        """
        v = min(0, self.left_pot)
        return v

    def max(self):
        """ Maximum value of the potential barrier for all positions and times.

            Returns:
                v: float
                    Maximum barrier height in reduced units
        """
        v = max(self.initial_height, self.final_height)
        return v

    def left_edge(self, t, E):
        """ Left edge of the barrier

            Position where a particle with the specified kinetic energy, 
            travelling from left to right, would collide with the barrier.

            Args:
                t: float
                    Time in reduced units
                E: float
                    Kinetic energy in reduced units

            Returns:
                x: float
                   Position of left edge in reduced units
        """
        x = 0
        return x

    def right_edge(self, t, E):
        """ Right edge of the barrier

            Position where a particle with the specified kinetic energy, 
            travelling from right to left, would collide with the barrier.

            Args:
                t: float
                    Time in reduced units
                E: float
                    Kinetic energy in reduced units

            Returns:
                x: float
                   Position of right edge in reduced units
        """
        x = self.width
        return x

    def _height(self, t):
        """ Height of the square barrier at the specified time.

            Args:
                t: float
                    Time in reduced units

            Returns:
                h: float
                    Barrier height in reduced units
        """
        if t < self.transition_time:
            h = self.initial_height
        else:
            h = self.final_height

        return h


class Square(SquareStep):
    """ Square barrier with option to suddenly change the barrier height.

        Args:
            params: dict
                Barrier parameters:
                    * width: width in fm
                    * initial_height: initial height in MeV  
                    * final_height: final height in MeV 
                    * transition_pos: wave packet classical position in fm when barrier changes height
                The first two are mandatory, the last two are optional.
            x_start: float 
                Wave packet start position in fm
            velocity: float
                Wave packet velocity in fm/s 

        Attributes:
            initial_height: float
                Initial height of the barrier in reduced units
            width: float
                Width of the barrier in reduced units
            transition_time: float
                Time at which barrier height changes
            final_height: float
                Final height of the barrier in reduced units
    """
    def __init__(self, params, x_start, velocity):
        super().__init__(params, x_start, velocity)


class BinnedBarrier():

    def __init__(self, barrier, bins, domain_size):

        self.barrier = barrier

        self.array = None
        self.prev_time = None

        dx = domain_size / bins
        self.x = np.arange(bins, dtype=np.float)
        self.x *= dx
        self.x -= 0.5 * (domain_size - dx)
        
    def get_array(self, t=0):

        b = self.barrier
        if b.has_changed(self.prev_time, t) or self.array is None:
            self.array = b.eval(self.x, t)

        self.prev_time = t
        return self.array

    def left_bin(self, t, E):
        xl = self.barrier.left_edge(t, E)
        bl = (xl - self.x[0]) / (self.x[1] - self.x[0])
        bl = int(bl)
        return bl

    def right_bin(self, t, E):
        xr = self.barrier.right_edge(t, E)
        br = (xr - self.x[0]) / (self.x[1] - self.x[0])
        br = int(br)
        return br

    def has_changed(self, t1, t2):
        return self.barrier.has_changed(t1, t2)

    def min(self):
        return self.barrier.min()

    def max(self):
        return self.barrier.max()


