import logging


class Logger:
    def __init__(self, log_file="all.log"):
        self.logger = self._setup_logger(log_file=log_file)

    def _setup_logger(self, log_file):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return logger
