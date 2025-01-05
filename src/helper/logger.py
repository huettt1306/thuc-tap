import logging
import os

def setup_logger(log_file="application.log", log_level=logging.INFO):
    """
    Thiết lập logger để ghi log vào một file.
    
    :param log_file: Đường dẫn file log.
    :param log_level: Mức độ log (mặc định: INFO).
    :return: Logger đã được cấu hình.
    """
    # Tạo thư mục cho log file nếu cần
    log_dir = os.path.dirname(log_file)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # Tạo logger
    logger = logging.getLogger(log_file)
    logger.setLevel(log_level)

    # Kiểm tra nếu logger đã được cấu hình trước đó
    if not logger.handlers:
        # Tạo handler để ghi log ra file
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)

        # Định dạng log
        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(formatter)

        # Thêm handler vào logger
        logger.addHandler(file_handler)

    return logger
