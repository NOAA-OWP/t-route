"""Basic Model Interface implementation for t-route."""

import numpy as np
import pandas as pd
from bmipy import Bmi
from pathlib import Path
import yaml

import nwm_routing.__main__ as tr


class bmi_troute(Bmi):

    def __init__(self):
        """Create a Bmi troute model that is ready for initialization."""
        super(bmi_troute, self).__init__()
        #self._model = None
        self._values = {}
        #self._var_units = {}
        self._var_loc = "node"
        self._var_grid_id = 0
        #self._grids = {}
        #self._grid_type = {}

        self._start_time = 0.0
        self._end_time = np.finfo("d").max
        self._time_units = "s"

    #----------------------------------------------
    # Required, static attributes of the model
    #----------------------------------------------
    _att_map = {
        'model_name':         'T-Route for Next Generation NWM',
        'version':            '',
        'author_name':        '',
        'grid_type':          'scalar', 
        'time_step_size':      1,       
        #'time_step_type':     'donno', #unused  
        #'step_method':        'none',  #unused
        #'time_units':         '1 hour' #NJF Have to drop the 1 for NGEN to recognize the unit
        'time_units':         'seconds' }

    #---------------------------------------------
    # Input variable names (CSDMS standard names)
    #---------------------------------------------
    _input_var_names = ['land_surface_water_source__volume_flow_rate',
                        'coastal_boundary__depth', #FIXME: this variable isn't a standard CSDMS name...couldn't find one more appropriate
                        'usgs_gage_observation__volume_flow_rate', #FIXME: this variable isn't a standard CSDMS name...couldn't find one more appropriate
                        'usace_gage_observation__volume_flow_rate', #FIXME: this variable isn't a standard CSDMS name...couldn't find one more appropriate
                        'rfc_gage_observation__volume_flow_rate' #FIXME: this variable isn't a standard CSDMS name...couldn't find one more appropriate
                       ]

    #---------------------------------------------
    # Output variable names (CSDMS standard names)
    #---------------------------------------------
    _output_var_names = ['channel_exit_water_x-section__volume_flow_rate',
                         'channel_water_flow__speed',
                         'channel_water__mean_depth',
                         'lake_water~incoming__volume_flow_rate',
                         'lake_water~outgoing__volume_flow_rate',
                         'lake_surface__elevation' #FIXME: this variable isn't a standard CSDMS name...couldn't find one more appropriate
                        ]

    #------------------------------------------------------
    # Create a Python dictionary that maps CSDMS Standard
    # Names to the model's internal variable names.
    #------------------------------------------------------
    _var_name_units_map = {
        'channel_exit_water_x-section__volume_flow_rate':['streamflow_cms','m3 s-1'],
        'channel_water_flow__speed':['streamflow_ms','m s-1'],
        'channel_water__mean_depth':['streamflow_m','m'],
        'lake_water~incoming__volume_flow_rate':['waterbody_cms','m3 s-1'],
        'lake_water~outgoing__volume_flow_rate':['waterbody_cms','m3 s-1'],
        'lake_surface__elevation':['waterbody_m','m'],
        #--------------   Dynamic inputs --------------------------------
        'land_surface_water_source__volume_flow_rate':['streamflow_cms','m3 s-1'],
        'coastal_boundary__depth':['depth_m', 'm'],
        'usgs_gage_observation__volume_flow_rate':['streamflow_cms','m3 s-1'],
        'usace_gage_observation__volume_flow_rate':['streamflow_cms','m3 s-1'],
        'rfc_gage_observation__volume_flow_rate':['streamflow_cms','m3 s-1']
    }

    #------------------------------------------------------
    # A list of static attributes. Not all these need to be used.
    #------------------------------------------------------
    _static_attributes_list = []


    #------------------------------------------------------------
    #------------------------------------------------------------
    # BMI: Model Control Functions
    #------------------------------------------------------------ 
    #------------------------------------------------------------

    #-------------------------------------------------------------------
    def initialize(self, bmi_cfg_file=None):

        args = tr._handle_args_v03(['-f', bmi_cfg_file])
        
        # -------------- Read in the BMI configuration -------------------------#
        bmi_cfg_file = Path(bmi_cfg_file)
        # ----- Create some lookup tabels from the long variable names --------#
        self._var_name_map_long_first = {long_name:self._var_name_units_map[long_name][0] for \
                                         long_name in self._var_name_units_map.keys()}
        self._var_name_map_short_first = {self._var_name_units_map[long_name][0]:long_name for \
                                          long_name in self._var_name_units_map.keys()}
        self._var_units_map = {long_name:self._var_name_units_map[long_name][1] for \
                                          long_name in self._var_name_units_map.keys()}
        
        # This will direct all the next moves.
        if bmi_cfg_file is not None:

            with bmi_cfg_file.open('r') as fp:
                cfg = yaml.safe_load(fp)
            self._cfg_bmi = self._parse_config(cfg)
        else:
            print("Error: No configuration provided, nothing to do...")
        
        # ------------- Initialize t-route model ------------------------------#
        (
            self._log_parameters, 
            self._preprocessing_parameters, 
            self._supernetwork_parameters, 
            self._waterbody_parameters, 
            self._compute_parameters, 
            self._forcing_parameters, 
            self._restart_parameters, 
            self._hybrid_parameters, 
            self._output_parameters, 
            self._parity_parameters, 
            self._data_assimilation_parameters,
        ) = tr._input_handler_v03(args)

        self._run_parameters = {
            'dt': self._forcing_parameters.get('dt'),
            'nts': self._forcing_parameters.get('nts'),
            'cpu_pool': self._compute_parameters.get('cpu_pool')
            }

        # Set number of time steps (1 hour)
        self._nts = 12
        
        # -------------- Initalize all the variables --------------------------# 
        # -------------- so that they'll be picked up with the get functions --#
        for var_name in list(self._var_name_units_map.keys()):
            # ---------- All the variables are single values ------------------#
            # ---------- so just set to zero for now.        ------------------#
            self._values[var_name] = np.zeros(self._network.dataframe.shape[0])
            setattr( self, var_name, 0 )
            
        '''
        # -------------- Update dimensions of DA variables --------------------# 
        ####################################
        # Maximum lookback hours from reservoir configurations
        usgs_shape = self._network._waterbody_types_df[self._network._waterbody_types_df['reservoir_type']==2].shape[0]
        usace_shape = self._network._waterbody_types_df[self._network._waterbody_types_df['reservoir_type']==3].shape[0]
        rfc_shape = self._network._waterbody_types_df[self._network._waterbody_types_df['reservoir_type']==4].shape[0]
        
        max_lookback_hrs = max(self._data_assimilation_parameters.get('timeslice_lookback_hours'),
                               self._waterbody_parameters.get('rfc').get('reservoir_rfc_forecasts_lookback_hours'))
        
        self._values['usgs_gage_observation__volume_flow_rate'] = np.zeros((usgs_shape,max_lookback_hrs*4))
        setattr( self, 'usgs_gage_observation__volume_flow_rate', 0 )
        self._values['usace_gage_observation__volume_flow_rate'] = np.zeros((usace_shape,max_lookback_hrs*4))
        setattr( self, 'usace_gage_observation__volume_flow_rate', 0 )
        self._values['rfc_gage_observation__volume_flow_rate'] = np.zeros((rfc_shape,max_lookback_hrs*4))
        setattr( self, 'rfc_gage_observation__volume_flow_rate', 0 )
        '''
        
        self._start_time = 0.0
        self._end_time = self._forcing_parameters.get('dt') * self._forcing_parameters.get('nts')
        self._time = 0.0
        self._time_step = self._forcing_parameters.get('dt')
        self._time_units = 's'
    
    def update(self):
        """Advance model by one time step."""
                
        # Set input data into t-route objects
        self._network._qlateral = pd.DataFrame(self._values['land_surface_water_source__volume_flow_rate'],
                                                            index=self._network.dataframe.index.to_numpy())
        self._network._coastal_boundary_depth_df = pd.DataFrame(self._values['coastal_boundary__depth'])


        # Create data assimilation object from da_sets for first loop iteration
        #TODO I'm not sure if this is the best way to do this. How will all of these variables be 
        # fed to t-route BMI?
        self._data_assimilation._usgs_df = pd.DataFrame(self._values['usgs_gage_observation__volume_flow_rate'])
        self._data_assimilation._last_obs_df = pd.DataFrame(self._values['lastobs__volume_flow_rate'])
        self._data_assimilation._reservoir_usgs_df = pd.DataFrame(self._values['reservoir_usgs_gage_observation__volume_flow_rate'])
        self._data_assimilation._reservoir_usgs_param_df = pd.DataFrame(self._values['reservoir_usgs__parameters'])
        self._data_assimilation._reservoir_usace_df = pd.DataFrame(self._values['reservoir_usace_gage_observation__volume_flow_rate'])
        self._data_assimilation._reservoir_usace_param_df = pd.DataFrame(self._values['reservoir_usace__parameters'])
        self._data_assimilation._da_parameter_dict = pd.DataFrame(self._values['DA_parameters'])

        ###NOTE: this is just a place holder, setting DA variables will be done with set_values...
        self._data_assimilation = tr.build_data_assimilation(
            self._network,
            self._data_assimilation_parameters,
            self._waterbody_parameters,
            [], #self._da_sets[0],
            self._forcing_parameters,
            self._compute_parameters)


        # Run routing
        (
            self._run_results, 
            self._subnetwork_list
        ) = tr.nwm_route(self._network.connections, 
                         self._network.reverse_network, 
                         self._network.waterbody_connections, 
                         self._network._reaches_by_tw,
                         self._compute_parameters.get('parallel_compute_method','serial'), 
                         self._compute_parameters.get('compute_kernel'),
                         self._compute_parameters.get('subnetwork_target_size'),
                         self._compute_parameters.get('cpu_pool'),
                         self._network.t0,
                         self._time_step,
                         self._nts,
                         self._forcing_parameters.get('qts_subdivisions', 12),
                         self._network.independent_networks, 
                         self._network.dataframe,
                         self._network.q0,
                         self._network._qlateral,
                         self._data_assimilation.usgs_df,
                         self._data_assimilation.lastobs_df,
                         self._data_assimilation.reservoir_usgs_df,
                         self._data_assimilation.reservoir_usgs_param_df,
                         self._data_assimilation.reservoir_usace_df,
                         self._data_assimilation.reservoir_usace_param_df,
                         self._data_assimilation.assimilation_parameters,
                         self._compute_parameters.get('assume_short_ts', False),
                         self._compute_parameters.get('return_courant', False),
                         self._network._waterbody_df,
                         self._waterbody_parameters,
                         self._network._waterbody_types_df,
                         self._network.waterbody_type_specified,
                         self._network.diffusive_network_data,
                         self._network.topobathy_df,
                         self._network.refactored_diffusive_domain,
                         self._network.refactored_reaches,
                         [], #subnetwork_list,
                         self._network.coastal_boundary_depth_df,
                         self._network.unrefactored_topobathy_df,)
        
        # update initial conditions with results output
        self._network.new_nhd_q0(self._run_results)
        self._network.update_waterbody_water_elevation()               
        
        # update t0
        self._network.new_t0(self._time_step,self._nts)

        # get reservoir DA initial parameters for next loop iteration
        self._data_assimilation.update(self._run_results,
                                       self._data_assimilation_parameters,
                                       self._run_parameters,
                                       self._network,
                                       [], #self._da_sets[run_set_iterator + 1]
                                      )
        
        (self._values['channel_exit_water_x-section__volume_flow_rate'], 
         self._values['channel_water_flow__speed'], 
         self._values['channel_water__mean_depth'], 
         self._values['lake_water~incoming__volume_flow_rate'], 
         self._values['lake_water~outgoing__volume_flow_rate'], 
         self._values['lake_surface__elevation'],
        ) = tr.create_output_dataframes(
            self._run_results, 
            self._nts, 
            self._network._waterbody_df,
            self._network.link_lake_crosswalk)
        
        self._time += self._time_step * self._nts

    def update_frac(self, time_frac):
        """Update model by a fraction of a time step.
        Parameters
        ----------
        time_frac : float
            Fraction fo a time step.
        """
        time_step = self.get_time_step()
        self._model.time_step = time_frac * time_step
        self.update()
        self._model.time_step = time_step

    def update_until(self, then):
        """Update model until a particular time.
        Parameters
        ----------
        then : float
            Time to run model until in seconds.
        """
        n_steps = (then - self.get_current_time()) / self.get_time_step()
        
        full_nts = self._nts
        self._nts = n_steps
        
        self.update()

        self._nts = full_nts - n_steps
        
        '''
        for _ in range(int(n_steps)):
            self.update()
        self.update_frac(n_steps - int(n_steps))
        '''

    def finalize(self):
        """Finalize model."""

        self._values = None
        self._var_loc = None
        self._var_grid_id = None
        self._var_name_map_long_first = None
        self._var_name_map_short_first = None
        self._var_units_map = None
        self._cfg_bmi = None
        self._network = None
        self._log_parameters = None
        self._preprocessing_parameters = None
        self._supernetwork_parameters = None
        self._waterbody_parameters = None
        self._compute_parameters = None
        self._forcing_parameters = None
        self._restart_parameters = None
        self._hybrid_parameters = None
        self._output_parameters = None
        self._parity_parameters = None
        self._data_assimilation_parameters = None
        self._run_parameters = None
        self._nts = None
        self._values = None
        self._start_time = None
        self._end_time = None
        self._time = None
        self._time_step = None
        self._time_units = None
        self._data_assimilation = None
        self._run_results = None
        self._subnetwork_list = None

    def get_var_type(self, var_name):
        """Data type of variable.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        str
            Data type.
        """
        return str(self.get_value_ptr(var_name).dtype)

    def get_var_units(self, var_name):
        """Get units of variable.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        str
            Variable units.
        """
        return self._var_units[var_name]

    def get_var_nbytes(self, var_name):
        """Get units of variable.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        int
            Size of data array in bytes.
        """
        return self.get_value_ptr(var_name).nbytes

    def get_var_itemsize(self, name):
        return np.dtype(self.get_var_type(name)).itemsize

    def get_var_location(self, name):
        return self._var_loc[name]

    def get_var_grid(self, var_name):
        """Grid id for a variable.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        int
            Grid id.
        """
        for grid_id, var_name_list in self._grids.items():
            if var_name in var_name_list:
                return grid_id

    def get_grid_rank(self, grid_id):
        """Rank of grid.
        Parameters
        ----------
        grid_id : int
            Identifier of a grid.
        Returns
        -------
        int
            Rank of grid.
        """
        return len(self._model.shape)

    def get_grid_size(self, grid_id):
        """Size of grid.
        Parameters
        ----------
        grid_id : int
            Identifier of a grid.
        Returns
        -------
        int
            Size of grid.
        """
        return int(np.prod(self._model.shape))

    def get_value_ptr(self, var_name):
        """Reference to values.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        array_like
            Value array.
        """
        return self._values[var_name]

    def get_value(self, var_name):
        """Copy of values.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        output_df : pd.DataFrame
            Copy of values.
        """
        output_df = self.get_value_ptr(var_name)
        return output_df

    def get_value_at_indices(self, var_name, dest, indices):
        """Get values at particular indices.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        dest : ndarray
            A numpy array into which to place the values.
        indices : array_like
            Array of indices.
        Returns
        -------
        array_like
            Values at indices.
        """
        dest[:] = self.get_value_ptr(var_name).take(indices)
        return dest

    def set_value(self, var_name, src):
        """
        Set model values
        
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        src : array_like
            Array of new values.
        """
        val = self.get_value_ptr(var_name)
        val[:] = src.reshape(val.shape)
        
        #self._values[var_name] = src

    def set_value_at_indices(self, name, inds, src):
        """Set model values at particular indices.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        src : array_like
            Array of new values.
        indices : array_like
            Array of indices.
        """
        val = self.get_value_ptr(name)
        val.flat[inds] = src

    def get_component_name(self):
        """Name of the component."""
        return self._name

    def get_input_item_count(self):
        """Get names of input variables."""
        return len(self._input_var_names)

    def get_output_item_count(self):
        """Get names of output variables."""
        return len(self._output_var_names)

    def get_input_var_names(self):
        """Get names of input variables."""
        return self._input_var_names

    def get_output_var_names(self):
        """Get names of output variables."""
        return self._output_var_names

    def get_grid_shape(self, grid_id, shape):
        """Number of rows and columns of uniform rectilinear grid."""
        var_name = self._grids[grid_id][0]
        shape[:] = self.get_value_ptr(var_name).shape
        return shape

    def get_grid_spacing(self, grid_id, spacing):
        """Spacing of rows and columns of uniform rectilinear grid."""
        spacing[:] = self._model.spacing
        return spacing

    def get_grid_origin(self, grid_id, origin):
        """Origin of uniform rectilinear grid."""
        origin[:] = self._model.origin
        return origin

    def get_grid_type(self, grid_id):
        """Type of grid."""
        return self._grid_type[grid_id]

    def get_start_time(self):
        """Start time of model."""
        return self._start_time

    def get_end_time(self):
        """End time of model."""
        return self._end_time

    def get_current_time(self):
        return self._time

    def get_time_step(self):
        return self._time_step

    def get_time_units(self):
        return self._time_units

    def get_grid_edge_count(self, grid):
        raise NotImplementedError("get_grid_edge_count")

    def get_grid_edge_nodes(self, grid, edge_nodes):
        raise NotImplementedError("get_grid_edge_nodes")

    def get_grid_face_count(self, grid):
        raise NotImplementedError("get_grid_face_count")

    def get_grid_face_nodes(self, grid, face_nodes):
        raise NotImplementedError("get_grid_face_nodes")

    def get_grid_node_count(self, grid):
        """Number of grid nodes.
        Parameters
        ----------
        grid : int
            Identifier of a grid.
        Returns
        -------
        int
            Size of grid.
        """
        return self.get_grid_size(grid)

    def get_grid_nodes_per_face(self, grid, nodes_per_face):
        raise NotImplementedError("get_grid_nodes_per_face")

    def get_grid_face_edges(self, grid, face_edges):
        raise NotImplementedError("get_grid_face_edges")

    def get_grid_x(self, grid, x):
        raise NotImplementedError("get_grid_x")

    def get_grid_y(self, grid, y):
        raise NotImplementedError("get_grid_y")

    def get_grid_z(self, grid, z):
        raise NotImplementedError("get_grid_z")
        
    def _parse_config(self, cfg):
        cfg_list = [cfg.get('flag'),cfg.get('file')]
        return cfg_list