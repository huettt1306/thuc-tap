import json

def get_trios_from_file(input_path):
    trios = {}
    with open(input_path, 'r') as f:
        # Bỏ qua dòng tiêu đề
        headers = f.readline()

        for line in f:
            line = line.strip()
            if not line:
                continue

            # Tách các trường thông tin từ dòng dữ liệu
            fields = line.split()
            if len(fields) < 8:
                continue

            family_id = fields[0]  # Family ID
            individual_id = fields[1]
            paternal_id = fields[2]
            maternal_id = fields[3]
            population = fields[6]  # Population

            # Chỉ lấy các bộ trio thuộc Population "KHV"
            if population == "KHV" and paternal_id != "0" and maternal_id != "0":
                if family_id not in trios:
                    trios[family_id] = {
                        "child": individual_id,
                        "mother": maternal_id,
                        "father": paternal_id
                    }

    return trios


# Đường dẫn file đầu vào và file kết quả đầu ra
input_path = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/integrated_call_samples_v3.20200731.ALL.ped"  # Thay bằng đường dẫn thật
output_path = "/home/huettt/Documents/nipt/NIPT-human-genetics/working/conf/trio.json"  # File JSON kết quả

# Lọc dữ liệu và lưu kết quả ra file JSON
trios = get_trios_from_file(input_path)
with open(output_path, 'w') as f:
    json.dump(trios, f, indent=4)

print(f"Trio data saved to {output_path}")
