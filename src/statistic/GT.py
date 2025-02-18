from helper.converter import convert_af_to_list
def get_af(gt, af):
    """
    Tính af cho biến thể trong kiểu gen gt
    0/0 -> 1-af
    0/1, 1/1 -> af
    Trả về -1 nếu kiểu gen không hợp lệ
    """
    allens = [int(a) for a in gt.split("/") if a != "."]
    if len(allens) != 2 or allens[0] < 0 or allens[1] > 1 or af < 0:
        return -1
    return int(100 * min(af, 1 - af))


def valid_gt(gt):
    """
    Kiểm tra 1 Kiểu gen có hợp lệ không
    """
    allens = [int(a) for a in gt.split("/") if a != "."]
    if len(allens) != 2:
        return False
    if allens[0] < 0 or allens[1] > 1:
        return False 
    return True


def has_gt(row, method):
    """
    Kiểm tra trong row, method có xác định ra kiểu gen hợp lệ không
    """ 
    return row[f"{method}"] and valid_gt(row[f"GT_{method}"])


def get_af_gt(row, method):
    """
    AF của kiểu gen mà method tìm được
    """
    if has_gt(row, method):
        return get_af(row[f"GT_{method}"], row["AF"])
    return -1


def get_af_gt_not_given(row, method, compare_with):
    """
    AF của Kiểu gen có trong method không có trong compare_with
    """
    if has_gt(row, method) and (not has_gt(row, compare_with)):
        return get_af(row[f"GT_{method}"], row["AF"])
    return -1

def get_af_gt_false(row, method, compare_with):
    """
    AF của kiểu gen có trong method khác với trong compare_with
    """
    if has_gt(row, method) and has_gt(row, compare_with) and row[f"GT_{method}"] != row[f"GT_{compare_with}"]:
        return get_af(row[f"GT_{method}"], row["AF"])
    return -1

def get_af_gt_true(row, method, compare_with):
    """
    AF của kiểu gen có trong method giống với trong compare_with
    """
    if has_gt(row, method) and has_gt(row, compare_with) and row[f"GT_{method}"] == row[f"GT_{compare_with}"]:
        return get_af(row[f"GT_{method}"], row["AF"])
    return -1


def get_af_gt_priv_true(row, method, truth, compare_with):
    """
    AF của kiểu gen method tìm được là đúng so với truth và khác so với compare_with
    """
    if has_gt(row, method) and row[f"GT_{method}"] == row[f"GT_{truth}"] and row[f"GT_{method}"] != row[f"GT_{compare_with}"]:
        return get_af(row[f"GT_{method}"], row["AF"])
    return -1

def get_af_gt_same_true(row, method, truth1, truth2):
    """
    AF của kiểu gen method tìm được là đúng so với truth1 và truth2
    """
    if has_gt(row, method) and row[f"GT_{method}"] == row[f"GT_{truth1}"] and row[f"GT_{method}"] == row[f"GT_{truth2}"]:
        return get_af(row[f"GT_{method}"], row["AF"])
    return -1

def get_af_gt_same_false(row, method, truth1, truth2):
    """
    AF của kiểu gen method tìm được là khác so với truth1 và truth2
    """
    if has_gt(row, method) and row[f"GT_{method}"] != row[f"GT_{truth1}"] and row[f"GT_{method}"] != row[f"GT_{truth2}"]:
        return get_af(row[f"GT_{method}"], row["AF"])
    return -1
