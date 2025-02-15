from helper.converter import convert_af_to_list
def valid_alt(gt):
    """
    Kiểm tra 1 Kiểu gen có alt không
    """
    allens = [int(a) for a in gt.split("/") if a != "."]
    if len(allens) != 2:
        return False
    return allens[1] == 1

def get_af(af):
    return int(100 * af)

def has_alt(row, method):
    """
    Kiểm tra trong row, method có xác định ra ALT không
    """ 
    return row[f"{method}"] and valid_alt(row[f"GT_{method}"])

def not_has_alt(row, method):
    """
    Kiểm tra trong row, method có xác định không có ALT không
    """ 
    return row[f"{method}"] and not valid_alt(row[f"GT_{method}"])


def get_af_alt(row, method):
    """
    AF của alt mà method tìm được
    """
    if has_alt(row, method):
        return get_af(row["AF"])
    return -1


def get_af_alt_not_given(row, method, compare_with):
    """
    AF của alt có trong method nhưng không có thông tin trong compare_with
    """
    if has_alt(row, method) and not row[f"{compare_with}"]:
        return get_af(row["AF"])
    return -1


def get_af_alt_false(row, method, compare_with):
    """
    AF của alt có trong method không có trong compare_with
    """
    if has_alt(row, method) and not_has_alt(row, compare_with):
        return get_af(row["AF"])
    return -1

def get_af_alt_true(row, method, compare_with):
    """
    AF của alt có trong method giống với trong compare_with
    """
    if has_alt(row, method) and has_alt(row, compare_with):
        return get_af(row["AF"])
    return -1


def get_af_alt_priv_true(row, method, truth, compare_with):
    """
    AF của alt method tìm được là đúng so với truth và khác so với compare_with
    """
    if has_alt(row, method) and has_alt(row, truth) and not_has_alt(row, compare_with):
        return get_af(row["AF"])
    return -1

def get_af_alt_same_true(row, method, truth1, truth2):
    """
    AF của alt method tìm được là đúng so với truth1 và truth2
    """
    if has_alt(row, method) and has_alt(row, truth1) and has_alt(row, truth2):
        return get_af(row["AF"])
    return -1

def get_af_alt_same_false(row, method, truth1, truth2):
    """
    AF của kiểu gen method tìm được là khác so với truth1 và truth2
    """
    if has_alt(row, method) and not_has_alt(row, truth1) and not_has_alt(row, truth2):
        return get_af(row["AF"])
    return -1
