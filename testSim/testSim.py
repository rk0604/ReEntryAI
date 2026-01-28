'''
Starting off with a 1D ODE -> x_dot = -x
# interval Euler: [x]_k+1=[x]_k + h⋅[f]([x]_k)
f(x) = -x
[f]([x]) = -[x]

Initial Uncertainty box: [x]_0 = [0.9, 1.1] --> simulating every number in this interval
step size (h) = 0.1
Next:
    1. [f]([x])
    2. multiply by h
    3. add to [x]_k

Write down:
'''
STEP_SIZE_H = 0.1
x_box = [0.9, 1.1] # [lo, hi]

# inclusion function version of f(x)
def inclusion_F(interval: list) -> list:
    lo = interval[0]
    hi = interval[1]
    return [(-1 * hi), (-1 * lo)] # [f]([x]) = -[x]

# interval addition helper
def interval_addition(interval1: list, interval2: list) -> list:
    # [a,b] + [c,d] = [a+c, b+d]
    lo_new = interval1[0] + interval2[0]
    hi_new = interval1[1] + interval2[1]
    return [lo_new, hi_new] # [a+c, b+d]

# scalar times interval helper
def scalar_times_interval(alpha: float, interval: list) -> list:
    # alpha * [lo, hi], flip bounds if alpha < 0
    lo = interval[0]
    hi = interval[1]
    if alpha >= 0:
        return [alpha * lo, alpha * hi]
    else:
        return [alpha * hi, alpha * lo]

results = {} # key = 1,2,3 and so on (1 = [x]_k). values are the intervals themselves

def main():
    global x_box  # we are updating the global x_box inside this function

    # store initial interval
    results[0] = x_box[:]
    print(f"[x]_0 : {x_box}")

    for i in range(10):
        # 1) [f]([x]_k)
        deriv_interval = inclusion_F(x_box)

        # 2) multiply by h
        update = scalar_times_interval(STEP_SIZE_H, deriv_interval)

        # 3) add to [x]_k
        new_step = interval_addition(x_box, update)

        results[i + 1] = new_step[:]
        print(f"[x]_{i+1} : {new_step}")

        x_box = new_step

if __name__ == "__main__":
    main()
'''
Sample Output:
[x]_0 : [0.9, 1.1]
[x]_1 : [0.79, 1.01]
[x]_2 : [0.6890000000000001, 0.931]
[x]_3 : [0.5959000000000001, 0.8621000000000001]
[x]_4 : [0.5096900000000001, 0.8025100000000001]
[x]_5 : [0.42943900000000007, 0.751541]
[x]_6 : [0.35428490000000007, 0.7085971]
[x]_7 : [0.28342519000000005, 0.67316861]
[x]_8 : [0.21610832900000004, 0.644826091]
[x]_9 : [0.15162571990000004, 0.6232152581]
[x]_10 : [0.08930419409000004, 0.60805268611]


'''
