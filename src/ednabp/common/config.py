from .default_settings import SETTINGS


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

    def __init__(
        self,
        logger=SETTINGS["config_basic"]["logger"],
        verbose=SETTINGS["config_basic"]["verbose"],
        dry=SETTINGS["config_basic"]["dry"],
    ):
        self.config_categories = ["basic"]
        self.logger = logger
        self.verbose = verbose
        self.dry = dry

    def __setattr__(self, name, value):
        if name == "verbose" and hasattr(self, "logger"):
            if value:
                self.logger.setLevel("INFO")
            else:
                self.logger.setLevel("WARNING")
        super().__setattr__(name, value)

    def get_basic_config(self) -> dict:
        """
        Retrieves the basic configuration settings.

        :returns: dict: A dictionary containing the verbose, dry run, and logger settings.
        """
        return {
            "logger": self.logger,
            "verbose": self.verbose,
            "dry": self.dry,
        }

    def add_machine_config(
        self,
        n_cpu=SETTINGS["config_machine"]["n_cpu"],
        memory=SETTINGS["config_machine"]["memory"],
    ):
        """
        Adds machine information to the configuration.

        :param n_cpu: Number of CPU cores allocated for the stage.
        :param memory: Amount of memory (in GB) allocated for the stage.
        """
        self.config_categories.append("machine")
        self.n_cpu = n_cpu
        self.memory = memory

    def get_machine_config(self) -> dict[str, int]:
        """
        Retrieves the machine information including CPU cores and memory.

        :returns: dict: A dictionary containing the number of CPU cores and amount of memory.
        """
        if "machine" not in self.config_categories:
            self.logger.warning(
                "Machine configuration not added. Run add_machine_config() first."
            )
            return
        return {"n_cpu": self.n_cpu, "memory": self.memory}

    # TODO(SW): Use these to help you organise all these parameters. dict[str:dict]
    def get_usearch_config(self) -> dict:
        return {}

    def get_denoise_config(self) -> dict:
        return {}

    def add_iqtree_config(self, model, boostrap, overwrite):
        self.config_categories.append("iqtree")
        self.iqtree_model = model
        self.iqtree_boostrap = boostrap
        self.iqtree_overwrite = overwrite

    def get_iqtree_config(self) -> dict:
        if "machine" not in self.config_categories:
            self.logger.warning(
                "Machine configuration not added. Run add_machine_config() first."
            )
            return
        if "iqtree" not in self.config_categories:
            self.logger.warning(
                "IQTree configuration not added. Run add_iqtree_config() first."
            )
            return
        return {
            "threads": self.n_cpu,
            "model": self.iqtree_model,
            "boostrap": self.iqtree_boostrap,
            "overwrite": self.iqtree_overwrite,
        }

    def add_cluster_config(self, reducer_kwargs, clusterer_kwargs, encode):
        self.config_categories.append("cluster")
        self.cluster_reducer_kwargs = reducer_kwargs
        self.cluster_clusterer_kwargs = clusterer_kwargs
        self.cluster_encode = encode

    def get_cluster_config(self) -> dict:
        if "cluster" not in self.config_categories:
            self.logger.warning(
                "Cluster configuration not added. Run add_cluster_config() first."
            )
            return
        return {
            "reducer_kwargs": self.cluster_reducer_kwargs,
            "clusterer_kwargs": self.cluster_clusterer_kwargs,
            "encode": self.cluster_encode,
        }

    def add_plot_config(self, show_plot, save_dir, overwrite):
        self.config_categories.append("plot")
        self.plot_show_plot = show_plot
        self.plot_save_dir = save_dir
        self.plot_overwrite = overwrite

    def get_plot_config(self) -> dict:
        if "plot" not in self.config_categories:
            self.logger.warning(
                "Plot configuration not added. Run add_plot_config() first."
            )
            return
        return {
            "show_plot": self.plot_show_plot,
            "save_dir": self.plot_save_dir,
            "overwrite": self.plot_overwrite,
        }
