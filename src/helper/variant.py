def convert_genotype(genotype):
    # Kiểm tra xem kiểu gen có đủ 2 phần tử và chúng là số nguyên
    if len(genotype) >= 2 and isinstance(genotype[0], int) and isinstance(genotype[1], int):
        # Đảm bảo sắp xếp theo thứ tự tăng dần
        sorted_genotype = sorted([genotype[0], genotype[1]])
        return f"{sorted_genotype[0]}/{sorted_genotype[1]}"
    return None

def get_af_variant(gt, af):
    """
    Tính af cho biến thể trong kiểu gen gt
    0/0 -> 1-af
    0/1, 1/1 -> af
    Trả về -1 nếu kiểu gen không hợp lệ
    """
    allens = [int(a) for a in gt.split("/") if a != "."]
    if len(allens) != 2 or allens[0] < 0 or allens[1] > 1:
        return -1
    if allens[1] == 0: 
        return 1 - af
    if allens[1] == 1:
        return af
    return -1

def check_valid_GT(gt):
    """
    Kiểm tra 1 Kiểu gen có hợp lệ không
    """
    allens = [int(a) for a in gt.split("/") if a != "."]
    if len(allens) != 2:
        return False
    if allens[0] < 0 or allens[1] > 1:
        return False 
    return True

def alt_variant(gt):
    """
    Kiểm tra 1 Kiểu gen có biến thể không
    """
    allens = [int(a) for a in gt.split("/") if a != "."]
    if len(allens) != 2:
        return False
    if allens[0] < 0 or allens[1] > 1:
        return False 
    return allens[1] == 1
    
def count_variant(row, method):
    """
    biến thể được liệt kê có trong cho mỗi method
    """
    if row[f"{method}"] and check_valid_GT(row[f"GT_{method}"]):
        return 1
    return 0 

def count_alt_variant(row, method, af):
    """
    biến thể hiếm mà method tìm được
    """
    if row[f"{method}"]:
        if get_af_variant(row[f"GT_{method}"], row["AF"]) >= af:
            return 1
    return 0 

def count_priv_gt(row, method, compare_with, af):
    """
    Đếm số lượng biến thể có trong method khác với biến thể trong compare_with
    """
    gt = row[f"GT_{method}"]
    if row[f"{method}"] and gt != row[f"GT_{compare_with}"] and get_af_variant(gt, row["AF"]) >= af:
        return 1
    return 0

def count_same_gt(row, method, compare_with, af):
    """
    Đếm số lượng biến thể có trong method giống với biến thể trong compare_with
    """
    gt = row[f"GT_{method}"]
    if row[f"{method}"] and gt == row[f"GT_{compare_with}"] and get_af_variant(gt, row["AF"]) >= af:
        return 1
    return 0

def count_priv_true_gt(row, method, truth, compare_with, af):
    """
    Biến thể method tìm được là đúng so với truth và khác so với compare_with
    """
    gt = row[f"GT_{method}"]
    if row[f"{method}"] and gt == row[f"GT_{truth}"] and gt != row[f"GT_{compare_with}"] and get_af_variant(gt, row["AF"]) >= af:
        return 1
    return 0

def count_same_true_gt(row, method, truth1, truth2, af):
    """
    Biến thể method tìm được là đúng so với truth1 và giống so với truth2
    """
    gt = row[f"GT_{method}"]
    if row[f"{method}"] and gt == row[f"GT_{truth1}"] and gt == row[f"GT_{truth2}"] and get_af_variant(gt, row["AF"]) >= af:
        return 1
    return 0

def count_same_false_gt(row, method, truth1, truth2, af):
    """
    Biến thể method tìm được là đúng so với truth1 và giống so với truth2
    """
    gt = row[f"GT_{method}"]
    if row[f"{method}"] and gt != row[f"GT_{truth1}"] and gt != row[f"GT_{truth2}"] and get_af_variant(gt, row["AF"]) >= af:
        return 1
    return 0
