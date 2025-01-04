import json

# Biến toàn cục
PATHS = {}
FTP_LINKS = {}
TRIO_DATA = {}
PARAMETERS = {}
TOOLS = {}

# Hàm load cấu hình
def load_config():
    global PATHS, FTP_LINKS, TRIO_DATA, PARAMETERS, TOOLS

    # Đọc file path.json
    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/path.json", "r") as file:
        PATHS = json.load(file)

    # Đọc file ftp.json
    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/ftp.json", "r") as file:
        FTP_LINKS = json.load(file)

    # Đọc file trio.json
    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/trio.json", "r") as file:
        TRIO_DATA = json.load(file)

    # Đọc file parameter.json
    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/parameter.json", "r") as file:
        PARAMETERS = json.load(file)

    # Đọc file tool_config.json
    with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/tool.json", "r") as file:
        TOOLS = json.load(file)

# Gọi hàm load_config() khi module được import
load_config()
