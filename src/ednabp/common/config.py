class Config:
    """
    Configuration class for stage settings.

    Attributes:
        verbose (bool): Indicates if verbose logging is enabled.
        dry (bool): Indicates if the stage should perform a dry run.
        logger (Logger): Logger object for logging messages.
        n_cpu (int): Number of CPU cores allocated for the stage.
        memory (int): Amount of memory (in GB) allocated for the stage.
    """

    def __init__(self, settings):
        self.config_categories = ["machine", "basic"]
        self.verbose = settings["verbose"]
        self.dry = settings["dry"]
        self.logger = settings["logger"]
        self.n_cpu = settings["n_cpu"]
        self.memory = settings["memory"]
        if self.verbose:
            self.logger.setLevel("INFO")
        else:
            self.logger.setLevel("WARNING")

    def __setattr__(self, name, value):
        if name == "verbose":
            if value:
                self.logger.setLevel("INFO")
            else:
                self.logger.setLevel("WARNING")
        super().__setattr__(name, value)

    def get_machine_info(self) -> dict[str, int]:
        """
        Retrieves the machine information including CPU cores and memory.

        :returns: dict: A dictionary containing the number of CPU cores and amount of memory.
        """
        return {"n_cpu": self.n_cpu, "memory": self.memory}

    def get_basic_configuration(self) -> dict:
        """
        Retrieves the basic configuration settings.

        :returns: dict: A dictionary containing the verbose, dry run, and logger settings.
        """
        return {
            "verbose": self.verbose,
            "dry_run": self.dry,
            "logger": self.logger,
        }

    # TODO(SW): Use these to help you organise all these parameters. dict[str:dict]
    def get_usearch_configuration(self) -> dict:
        return {}

    def get_denoise_configuration(self) -> dict:
        return {}

    def add_iqtree_configuration(self, model, boostrap, overwrite):
        self.config_categories.append("iqtree")
        self.iqtree_model = model
        self.iqtree_boostrap = boostrap
        self.iqtree_overwrite = overwrite

    def get_iqtree_configuration(self) -> dict:
        if "iqtree" not in self.config_categories:
            self.logger.warning(
                "WARNING: IQTree configuration not added. Run add_iqtree_configuration() first."
            )
            return
        return {
            "threads": self.n_cpu,
            "model": self.iqtree_model,
            "boostrap": self.iqtree_boostrap,
            "overwrite": self.iqtree_overwrite,
        }

    def add_seqcluster_configuration(
        self, reducer_kwargs, clusterer_kwargs, encode
    ):
        self.config_categories.append("seqcluster")
        self.seqclu_reducer_kwargs = reducer_kwargs
        self.seqclu_clusterer_kwargs = clusterer_kwargs
        self.seqclu_encode = encode

    def get_seqcluster_configuration(self) -> dict:
        if "seqcluster" not in self.config_categories:
            self.logger.warning(
                "WARNING: SeqCluster configuration not added. Run add_seqcluster_configuration() first."
            )
            return
        return {
            "reducer_kwargs": self.seqclu_reducer_kwargs,
            "clusterer_kwargs": self.seqclu_clusterer_kwargs,
            "encode": self.seqclu_encode,
        }
