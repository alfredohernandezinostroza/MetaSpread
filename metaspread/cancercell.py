import mesa
from metaspread.vessel import Vessel

class CancerCell(mesa.Agent):

    def __init__(self, unique_id, model, grid, grid_id, phenotype, ecm, mmp2):
        super().__init__(unique_id, model)
        self.grid = grid
        self.grid_id = grid_id
        self._apply_phenotype(phenotype)
        self.ecm = ecm
        self.mmp2 = mmp2
        self.agent_type = "cell"
        self.ruptured = False #need to be able do use data collector on agents

    def _apply_phenotype(self, phenotype):
        """Set the phenotype and its associated motility parameters.

        Shared by __init__ and by EMT/MET switching so the diffusion coefficient
        and haptotaxis sensitivity always match the current phenotype.
        """
        self.phenotype = phenotype
        if phenotype == "mesenchymal":
            self.diff_coeff = self.model.config.dM
            self.phi = self.model.config.phiM
        else:
            self.diff_coeff = self.model.config.dE
            self.phi = self.model.config.phiE

    def step(self): #what will the agent do every time a step is made
        self.update_phenotype()
        self.move()

    def update_phenotype(self):
        """Optionally switch phenotype (EMT/MET plasticity).

        Off by default: with emt_prob=met_prob=0 and enable_hypoxia_emt=False the
        method returns without consuming the RNG, so default runs are unchanged.
        - stochastic: epithelial->mesenchymal (EMT) with prob emt_prob,
          mesenchymal->epithelial (MET) with prob met_prob, every step.
        - hypoxia (enable_hypoxia_emt, requires enable_oxygen): the EMT drive only
          applies where local oxygen is below hypoxia_threshold.
        """
        cfg = self.model.config
        emt_p = cfg.emt_prob
        met_p = cfg.met_prob
        if cfg.enable_hypoxia_emt:
            x, y = self.pos
            oxygen_here = self.model.oxygen[self.grid_id - 1][0, x, y]
            if oxygen_here >= cfg.hypoxia_threshold:
                emt_p = 0.0  # sufficient oxygen: no hypoxia-driven EMT here
        if self.phenotype == "epithelial":
            if emt_p > 0 and self.random.random() < emt_p:
                self._apply_phenotype("mesenchymal")
        else:  # mesenchymal
            if met_p > 0 and self.random.random() < met_p:
                self._apply_phenotype("epithelial")

    def move(self):
        #fixed probabilities can be given to fix the movement of the cells towards a certain direction
        #if no fixed probabilities are given, the probabilities are calculated based on the ECM and MMP-2 concentration
        #(the fixed probabilities are only used for testing)
        fixed_p_left    = self.model.fixed_p_left
        fixed_p_right   = self.model.fixed_p_right
        fixed_p_top     = self.model.fixed_p_top
        fixed_p_bottom  = self.model.fixed_p_bottom
        th = self.model.config.th
        xh = self.model.config.xh
        time = self.model.schedule.time
        possible_steps = self.grid.get_neighborhood(
            self.pos,
            moore=False,
            include_center=True)
        x, y = self.pos
        on_left_border    = self.grid.out_of_bounds((x-1,y))
        on_right_border   = self.grid.out_of_bounds((x+1,y))
        on_top_border     = self.grid.out_of_bounds((x,y+1))
        on_bottom_border  = self.grid.out_of_bounds((x,y-1))
        p_left   = (fixed_p_left if fixed_p_left is not None
                        else 0 if on_left_border 
                        else (th/xh**2*(self.diff_coeff-self.phi/4*
                            (0  if on_right_border
                                else self.ecm[0,x+1,y]-self.ecm[0,x-1,y]))))
        p_right  = (fixed_p_right if fixed_p_right is not None
                        else 0 if on_right_border 
                        else (th/xh**2*(self.diff_coeff+self.phi/4*
                            (0  if on_left_border
                                else self.ecm[0,x+1,y]-self.ecm[0,x-1,y]))))
        p_top    = (fixed_p_top if fixed_p_top is not None
                        else 0 if on_top_border 
                        else (th/xh**2*(self.diff_coeff+self.phi/4*
                            (0  if on_bottom_border
                                else self.ecm[0,x,y+1]-self.ecm[0,x,y-1]))))
        p_bottom = (fixed_p_bottom if fixed_p_bottom is not None
                        else 0 if on_bottom_border 
                        else (th/xh**2*(self.diff_coeff-self.phi/4*
                            (0  if on_top_border
                                else self.ecm[0,x,y+1]-self.ecm[0,x,y-1]))))
        p_stay      = 1-(p_left+p_right+p_top+p_bottom)

        weights=[]
        for x2,y2 in possible_steps:
            if x2 < x:
                weights.append(p_left)
            elif x2>x:
                weights.append(p_right)
            elif y2<y:
                weights.append(p_bottom)
            elif y2>y:
                weights.append(p_top)
            else:
                weights.append(p_stay)
        #there is a chance that the probability should be calculated at the end as 1-sum(weights), have to check


        # new_position = (x,y+1)
        new_position = self.random.choices(possible_steps,weights,k=1)[0]
        is_vessel = False
        is_ruptured = False
        for agent in self.grid.get_cell_list_contents([(new_position)]):
            if isinstance(agent, Vessel):
                is_ruptured = agent.ruptured
                is_vessel=True
                break
        if is_vessel and self.grid_id == 1 and (is_ruptured or self.phenotype == "mesenchymal"): 
                x, y = new_position
                on_left_border    = self.grid.out_of_bounds((x-1,y))
                on_right_border   = self.grid.out_of_bounds((x+1,y))
                on_top_border     = self.grid.out_of_bounds((x,y+1))
                on_bottom_border  = self.grid.out_of_bounds((x,y-1))
                mesenchymal_ccells_to_travel  = [agent for agent in self.grid.get_cell_list_contents([(x,y)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_left_border   else [agent for agent in self.grid.get_cell_list_contents([(x-1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_right_border  else [agent for agent in self.grid.get_cell_list_contents([(x+1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_top_border    else [agent for agent in self.grid.get_cell_list_contents([(x,y-1)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_bottom_border else [agent for agent in self.grid.get_cell_list_contents([(x,y+1)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                epithelial_ccells_to_travel  = [agent for agent in self.grid.get_cell_list_contents([(x,y)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_left_border   else [agent for agent in self.grid.get_cell_list_contents([(x-1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_right_border  else [agent for agent in self.grid.get_cell_list_contents([(x+1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_top_border    else [agent for agent in self.grid.get_cell_list_contents([(x,y-1)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_bottom_border else [agent for agent in self.grid.get_cell_list_contents([(x,y+1)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
        
                #if there are not clusters at that time in the vasculature dict, create a new key for that time
                #and add the tuple

                vasculature_time = self.model.config.vasculature_time
                if self.model.vasculature.get(time + vasculature_time,False):
                    self.model.vasculature[time + vasculature_time] += [(len(mesenchymal_ccells_to_travel), len(epithelial_ccells_to_travel))]
                # if there are clusters, add the tuple to that key
                else:
                    self.model.vasculature[time + vasculature_time] = [(len(mesenchymal_ccells_to_travel), len(epithelial_ccells_to_travel))]
                for ccell in mesenchymal_ccells_to_travel + epithelial_ccells_to_travel:
                    ccell.grid.remove_agent(ccell)
                    ccell.model.schedule.remove(ccell)
        else:
            if self.model.config.carrying_capacity > len([cell for cell in self.grid.get_cell_list_contents([new_position]) if cell.agent_type == 'cell']):
                self.grid.move_agent(self, new_position)
