import mesa


class ImmuneCell(mesa.Agent):
    """A spatial immune (effector) cell performing surveillance and cytotoxicity.

    Phase 1 immune model: immune cells diffuse and kill nearby cancer cells with
    probability ``immune_kill_prob``. Systemic effector/suppressor recruitment
    dynamics (the digital-twin ODE layer) are deferred to a later phase. Only
    created when ``enable_immune`` is set, so default simulations are unaffected.

    Exposes ``agent_type``/``phenotype``/``ruptured`` so the model's DataCollector
    (which reports those fields for every agent) works unchanged.
    """

    def __init__(self, unique_id, model, grid, grid_id):
        super().__init__(unique_id, model)
        self.grid = grid
        self.grid_id = grid_id
        self.agent_type = "immune"
        self.phenotype = False  # needed for the data collector
        self.ruptured = False   # needed for the data collector

    def step(self):
        self.move()
        self.attack()

    def move(self):
        """Diffusive random walk (same th/xh scaling as cancer cell diffusion)."""
        cfg = self.model.config
        possible_steps = self.grid.get_neighborhood(
            self.pos, moore=False, include_center=True)
        p = cfg.th / cfg.xh ** 2 * cfg.immune_diff_coeff
        n_neighbors = len(possible_steps) - 1
        stay = 1 - n_neighbors * p
        if stay < 0:  # diffusion too high for a proper distribution: move uniformly
            weights = [1.0] * len(possible_steps)
        else:
            weights = [stay if step == self.pos else p for step in possible_steps]
        # Device geometry (Phase 2): never step into a wall cell.
        if self.model.wall_mask is not None:
            weights = [0.0 if self.model.wall_mask[tuple(step)] else w
                       for step, w in zip(possible_steps, weights)]
            if sum(weights) <= 0:  # walled in: stay put
                weights = [1.0 if step == self.pos else 0.0 for step in possible_steps]
        new_position = self.random.choices(possible_steps, weights, k=1)[0]
        self.grid.move_agent(self, new_position)

    def attack(self):
        """Kill cancer cells in the Moore neighborhood with immune_kill_prob."""
        kill_prob = self.model.config.immune_kill_prob
        if kill_prob <= 0:
            return
        neighborhood = self.grid.get_neighborhood(
            self.pos, moore=True, include_center=True)
        targets = [a for a in self.grid.get_cell_list_contents(neighborhood)
                   if a.agent_type == "cell"]
        for target in targets:
            if self.random.random() < kill_prob:
                self.grid.remove_agent(target)
                self.model.schedule.remove(target)
                self.model.cancer_cells_counter[self.grid_id - 1] -= 1
