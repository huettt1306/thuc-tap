import json
from ftplib import FTP

# Kết nối đến FTP server
ftp_server = 'ftp.sra.ebi.ac.uk'
ftp = FTP(ftp_server)
ftp.login()

# Đường dẫn đến thư mục gốc
base_path = '/vol1/run/'
result_path = '/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/test.txt'

# Đọc tên mẫu từ file JSON hoặc bất kỳ cấu trúc dữ liệu nào
with open("/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/trio.json", "r") as file:
        sample_data = json.load(file)

# Chuyển đổi tất cả các tên mẫu (cả child, mother, father) thành một set để kiểm tra nhanh hơn
sample_names = set()
for sample_group in sample_data.values():
    sample_names.update(sample_group.values())  # Lấy tất cả giá trị 'child', 'mother', 'father'

print (sample_names)

# Lưu tên mẫu đã có đường dẫn vào set
found_samples = set()

target_dirs = {
    #'ERR398': (3988970, 3989050),  # Duyệt từ ERR3988970 đến ERR3980950
    'ERR324': (3242356, 3242356),  # Duyệt từ ERR3242000 đến ERR3242050
}

# Hàm kiểm tra và lưu đường dẫn FTP vào file
def check_and_save_paths():
    with open(result_path, 'w') as f_out:
        # Duyệt qua thư mục cấp 1 (A)
        for dir_a, (start_number, end_number) in target_dirs.items():
            print(f"Đang duyệt thư mục chính: {dir_a}")
            ftp.cwd(base_path)
            ftp.cwd(dir_a)

            # Duyệt qua các thư mục con theo thứ tự tăng dần bắt đầu từ số đã chỉ định
            for current_number in range(start_number, end_number + 1):
                dir_b = f"ERR{current_number:07d}"  # Tạo tên thư mục theo định dạng ERRxxxxxxx
                print(dir_b)
                try:
                    ftp.cwd(dir_b)  # Thử vào thư mục
                    # Duyệt qua các file trong thư mục B
                    for filename in ftp.nlst():
                        if filename.endswith(".final.cram"):
                            # Tách phần đầu tên file (trước dấu '_')
                            parts = filename.split('.')
                            if parts:
                                sample_name = parts[0]  # Tên mẫu là phần đầu tiên trước dấu '_'
                                # Kiểm tra nếu tên mẫu hợp lệ và có trong dữ liệu JSON
                                print(sample_name)
                                if sample_name in sample_names:
                                    remote_path = f"{ftp_server}{base_path}{dir_a}/{dir_b}/{filename}"
                                    # Ghi đường dẫn FTP vào file
                                    f_out.write(f"{remote_path}\n")
                                    found_samples.add(sample_name)  # Thêm mẫu vào set đã tìm thấy
                            else:
                                print(f"Bỏ qua file không hợp lệ: {filename}")
                    ftp.cwd(base_path + dir_a)  # Quay lại thư mục chính
                except Exception as e:
                    # Nếu thư mục không tồn tại, dừng duyệt
                    print(f"Không tìm thấy thư mục {dir_b}, dừng duyệt.")
                    break

                current_number += 1  # Tăng số thứ tự thư mục con

    print("Đã lưu các đường dẫn vào ftppath.txt.")


# Gọi hàm để kiểm tra và lưu đường dẫn FTP
check_and_save_paths()

# Kiểm tra các mẫu trong JSON chưa có đường dẫn
missing_samples = sample_names - found_samples
if missing_samples:
    print("Các mẫu chưa có đường dẫn tải:")
    for missing in missing_samples:
        print(missing)
else:
    print("Tất cả các mẫu đều đã có đường dẫn tải.")

# Đóng kết nối FTP
ftp.quit()
