def calculate_scale(old_lower, old_upper, new_lower, new_upper):
    old_lower, old_upper, new_lower, new_upper = map(float, [old_lower, old_upper, new_lower, new_upper])
    delta_old = old_upper - old_lower
    delta_new = new_upper - new_lower
    scale = delta_new / delta_old
    return scale


def transform_dimension(value, translate, scale):
    return (value - translate) * scale


def transform_dimension_inverse(value, translate, scale):
    return (value / scale) + translate


def _test_one():
    old_lower, old_upper = -180, 180
    new_lower, new_upper = 0, 2**21
    scale = calculate_scale(old_lower, old_upper, new_lower, new_upper)
    translate = old_lower
    val = 0
    x = transform_dimension(val, translate, scale)
    y = transform_dimension_inverse(x, translate, scale)
    print val, '=>', x, '=>', y


def _test_brute():
    old_lower, old_upper = -180, 180
    new_lower, new_upper = 0, 2**21
    scale = calculate_scale(old_lower, old_upper, new_lower, new_upper)
    translate = old_lower
    # go both forward and backward
    val = 0
    for i in range(old_lower, old_upper + 1):
        x = transform_dimension(val, translate, scale)
        y = transform_dimension_inverse(x, translate, scale)
        assert y == val


if __name__ == "__main__":
    _test_one()
    _test_brute()


