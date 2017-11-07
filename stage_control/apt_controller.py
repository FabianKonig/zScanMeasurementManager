if __name__ == '__main__':
    import thorlabs_apt as apt

else:
    from . import thorlabs_apt as apt

import numpy as np
from math import isclose


# Check that both stages are connected to the PC and are detected
assert len(apt.list_available_devices()) == 2




class APT_Controller:
    """ Main class to control the stages. This class basically wraps functionality from
        the thorlabs_apt package in functions that are useful for our purposes.
    """

    def __init__(self):
        """ Store all available motors as apt.Motor instances, change their velocity parameters,
            set their maximum position, calibrate them by moving them to their home position and
            finally, move them to their maximum position. The last step is required in order to
            start the measurement in front of the beam waist of the focused beam (The home position
            of the stages is behind the focal spot and the maximum positions of the stages are in
            front of it.)
        """
        self._motors = self.get_motors()
        self.total_travel_distance = None  # The total distance that all stages combined can travel.
                                           # Is set in self.set_stages_maximum_positions()
        self.combined_position = self._get_combined_position()
                                           # The current position of both stages combined. 
                                           # Is updated after each movement.



        for motor in self._motors:
            assert motor.get_stage_axis_info()[2] == 1  # assert that unit is millimetres
            motor.set_velocity_parameters(0, 3.5, 2.1)
        
        self.set_stages_maximum_positions()
        
        if False:#self.combined_position > 0.1:  # move the motors quickly close to the home position.
            for motor in self._motors:
                motor.move_to(0.01)
            self.combined_position = self._get_combined_position
            self.block_while_moving_then_assert_and_update_position(0.02)

        #self.move_stages_home()
        self.initialise_stages()  # Currently, this just invokes move_stages_to_maximum_positions()


    def get_motors(self):
        """ Returns all USB connected motors as an array of apt.Motor instances """

        devices = apt.list_available_devices()
        return np.array([apt.Motor(devices[i][1]) for i in np.arange(len(devices))])


    def _get_combined_position(self):
        """ Returns the position of all stages combined """

        combined_position = 0
        for motor in self._motors:
            combined_position += motor.position

        return combined_position


    def set_stages_maximum_positions(self):
        """ Set the maximum positions of the stages below 25mm. I have empirically discovered
            that they cannot move as far as 25mm, unfortunately.
        """
        
        for motor in self._motors:
            # Store any other stage info:
            axis_info = motor.get_stage_axis_info()
            max_pos = axis_info[1]

            if motor.serial_number == 83822061:
                max_pos = 22
            elif motor.serial_number == 83825266:
                max_pos = 23

            motor.set_stage_axis_info(axis_info[0], max_pos, axis_info[2], axis_info[3])

        # Set the instance property self.total_travel_distance
        self.total_travel_distance = 0
        
        for motor in self._motors:
            self.total_travel_distance += motor.get_stage_axis_info()[1]


    def block_while_moving_then_assert_and_update_position(self, expected_new_combined_position):
        """ Function ensuring that the program halts until the stages are no longer in motion. 
            After movement has completed, the new combined position is checked against the 
            function parameter "expected_new_combined_position". If the check succeeds, 
            self.combined_position is updated.

        Input:
            expected_new_combined_position  The new combined stage position to be expected after
                                            stage movement is completed.
        """

        while(True):
            break_condition = True
            for motor in self._motors:
                break_condition *= not motor.is_in_motion
            if(break_condition):
                break

        new_combined_position = 0
        for motor in self._motors:
            new_combined_position += motor.position

        # Assert the position with a precision of up to 0.05mm = 50µm.
        assert isclose(new_combined_position, expected_new_combined_position, abs_tol=0.05)
        self.combined_position = new_combined_position


    def move_stages_home(self):
        for motor in self._motors:
            motor.move_home()

        self.block_while_moving_then_assert_and_update_position(0)


    def move_stages_to_maximum_positions(self):
        for motor in self._motors:
            max_position = motor.get_stage_axis_info()[1]
            motor.move_to(max_position)

        self.block_while_moving_then_assert_and_update_position(self.total_travel_distance)


    def initialise_stages(self):
        """ Currently, the stages' initial position is their maximum position. """
        self.move_stages_to_maximum_positions()


    def move_in_steps(self, tot_num_of_pos, direction):
        if direction == "forward":
            direction = +1
        elif direction == "backward":
            direction = -1
        else:
            raise Exception("Direction must be a string of either 'forward' or 'backward'.")

        expected_new_combined_position = 0

        for motor in self._motors:
            min_stage_position = motor.get_stage_axis_info()[0]  # Das als static wäre gut
            max_stage_position = motor.get_stage_axis_info()[1]

            # The first position is always the initial position. Hence "tot_num_of_pos-1"
            step_size = direction * max_stage_position / (tot_num_of_pos-1)
            new_position = motor.position + step_size

            # new_position is only allowed to exceed the minimum/maximum position by 10µm
            if new_position > max_stage_position:
                assert isclose(new_position, max_stage_position, abs_tol=0.01)
                new_position = max_stage_position
            elif new_position < min_stage_position:
                assert isclose(new_position, min_stage_position, abs_tol=0.01)
                new_position = min_stage_position

            expected_new_combined_position += new_position
            motor.move_to(new_position)

        self.block_while_moving_then_assert_and_update_position(expected_new_combined_position)



if __name__ == '__main__':
    aptc = APT_Controller()

    # Check if maximum stage positions have been properly set:
    total_max_stage_positions = 0
    for motor in aptc._motors:
        total_max_stage_positions += motor.get_stage_axis_info()[1]

    assert isclose(total_max_stage_positions, 22+23, abs_tol=0.05)

    total_steps = range(45)

    import time
    start = time.time() 
    for i in total_steps:
        aptc.move_in_steps(len(total_steps), "backward")
    print((time.time() - start) / 60)
