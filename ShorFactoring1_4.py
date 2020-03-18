from math import log2, e, pi, ceil
import random
import numpy as np
np.set_printoptions(threshold=np.inf)  

#a^b%c，快速模运算
def quick_mod(a, b, c):
    a = a % c
    bin_b = bin(b)[2:][::-1]
    bin_b_list = []
    bin_b_list.append(a)
    for _ in range(len(bin_b)-1):
        a = (a*a) % c
        bin_b_list.append(a)
    result = 1
    for i,flag in enumerate(bin_b):
        if flag == '1':
            result = (result*bin_b_list[i])%c
    return result

# 计算a^b次方，复杂度O(log(b))
def quick_exponential_operation(a,b):
    bin_b = bin(b)[2:][::-1]
    bin_b_list = []
    bin_b_list.append(a)
    for _ in range(len(bin_b)-1):
        a = (a*a) 
        bin_b_list.append(a)
    result = 1
    for i,flag in enumerate(bin_b):
        if flag == '1':
            result = (result*bin_b_list[i])
    return result

# Miller-Rabin算法判断质数
# num大于等于2
def is_prime(num):
    if num == 2:
        return True
    if num % 2 == 0:
        return False
    if num == 3:
        return True
    k = 0
    temp = num - 1
    while temp % 2 == 0:
        k+=1
        temp //= 2
    q = temp
    #10次测试，一个非素整数通过的概率为10^-6
    for _ in range(10):
        a = random.randint(2, num-2)
        if quick_mod(a, q, num) == 1:
            continue
        for j in range(k):
            if quick_mod(a, (q*2**j), num) == (num - 1):
                break
        else:
            #若执行到此，不满足素数的必要条件，返回False
            return False
    #通过十次检测，返回True
    return True

# N大于等于2
def is_prime_power(N):
    # N是k^p形式，k是质数，p大于等于1
    for i in range(1, int(log2(N))+1):
        # 如果能表示成某数k的p次方，那么当k不为质数时，N当然不会为某个质数d的次方
        if abs(int(N**(1/i)) - N**(1/i)) < 0.00001:
            if is_prime(int(N**(1/i))):
                return True, (int(N**(1/i)), i)
    return False, -1

# gcd(a,b) = gcd(b, a mod b)
def my_gcd(a, b):
    if a < b:
        a,b = b,a
    if b == 0:
        return a
    else:
        return my_gcd(b, a%b)

# 矩阵状态
class QuantumState():
    # 传入的状态列表，来实现初始化
    def __init__(self, form, parameter):
        # 可以用bit串形式对于的一个基向量来初始化
        if form == 'bit':
            temp = int(parameter, 2)
            self.state = np.zeros(2**len(parameter), dtype='complex64')
            self.state[temp] = 1
        # 只有在使用基本门操作生成新的量子态的时候才可以用'array'构造
        # array为numpy类型
        # 物理实现中，一般只能直接获得一些简单的量子态
        elif form == 'array':
            self.state = parameter
        # n个为0的qubit叠加形成的态
        elif form == 'zeros':
            self.state = np.zeros(2**parameter, dtype='complex64')
            self.state[0] = 1
        else:
            print('Error!')
    def size(self):
        return self.state.shape[0]
    # 观测
    # 观测qubit_list中的相应位置的qubit，然后状态坍缩
    # 先计算qubit_list中所有情况中每一个情况对应的概率，然后根据该概率进行坍缩
    # 坍缩的结果为，只取一种结果，即母串中只保留和有该结果（子串的）的情形
    def observe(self, qubit_list):
        if qubit_list == 'all':
            qubit_list = [i for i in range(int(log2(self.state.shape[0])))]
        if abs(((abs(self.state)**2).sum() - 1)) > 0.000001:
            print('Error!', '向量二范数不为1！')
            print((abs(self.state)**2).sum())
        length = self.state.shape[0]
        probability_list = [0]*(2**len(qubit_list))
        for i in range(length):
            bin_i = bin(i)[2:]
            # 补零
            bin_i = '0'*(int(log2(length))-len(bin_i))+bin_i
            # 抽取特定位置的bit
            extract_bit = [bin_i[i] for i in qubit_list]
            q = extract_bit[0]
            for q_ in extract_bit[1:]:
                q += q_
            extract_number = int(q, 2)
            # print(extract_number)
            # print((abs(self.state[i])**2))
            probability_list[extract_number] += (abs(self.state[i])**2)
        random_num = random.random()
        temp = 0
        count = 0
        # 计算会有误差，如果random_num为1，那么计算出来的temp也有可能小于1
        # 所以要对count做出判断，以防越界的情况发生
        while count < len(probability_list) and temp < random_num:
            temp += probability_list[count]
            count += 1
        count -= 1
        result_state = np.zeros(self.state.shape[0], dtype = np.complex_)
        # 此时，母串中若含count对应的bit串，那么则保留，否则一律归零
        # 然后进行归一化操作
        for i in range(length):
            bin_i = bin(i)[2:]
            # 补零
            bin_i = '0'*(int(log2(length))-len(bin_i))+bin_i
            # 抽取特定位置的bit
            extract_bit = [bin_i[i] for i in qubit_list]
            q = extract_bit[0]
            for q_ in extract_bit[1:]:
                q += q_
            extract_number = int(q, 2)
            if extract_number == count:
                result_state[i] = self.state[i]
            # print(abs(result_state).sum())
        result_state /= ((abs(result_state)**2).sum()**0.5)
        count = bin(count)[2:]
        # 补零
        count = '0'*(int(len(qubit_list))-len(count))+count
        print(count)
        # 返回坍缩后的状态以及观测的qubit串坍缩后形成的状态
        return QuantumState('array', result_state), QuantumState('bit', count)
    def kron(self, state):
        return QuantumState('array', np.kron(self.state, state))

    # 切片处理，即只取几位的qubit来研究
    # 注意，只提取的几位qubit，仍可能和未提取的qubit有纠缠
    # 对于m位qubit串，先提取前a位，再提取后m-a位，二者得到的结果进行张量积运算，不一定得到原qubit串
    # 如果二者不纠缠，那么便可以得到原qubit串 
    def front_slice(self):
        pass
    # 分解为最简的张量积形式
    # 返回QuantumState组成的列表，按迭代的顺序对其中的QuantumState进行张量积运算，就得到该state
    # 可惜的是，复数域的分解并不是唯一的
    # 比如，对于(a0|0>+a1|1>)张量乘(b0|0>+b1|1>)=a0b0|00>+a0b1|01>+a1b0|10>+a1b1|11>
    # 可以看到，若a0和a1两个复数，顺时针旋转theta度，而b0和b1两个复数，逆时针旋转theta度，那么等式右侧是不会变的
    # 而等式左侧变化
    # 也就是说，分解并不唯一。更具体来说，分解的幅度是确定的，但是旋转角便是不确定的了
    def decompose(self):
        pass 
    

    # 反张量积运算
    # 设本类为c，参数state为b，有a张量积b=c，本函数返回a
    def rev_kron(self, state):
        if abs(state.state).sum() == 0:
            print('Error')
            return 
        count = 0
        # 找到state中第一个不为零的元素
        while state.state[count] == 0:
            count += 1
        result_state = np.zeros(self.state.shape[0]//state.state.shape[0], dtype = np.complex_)
        for i in range(result_state.shape[0]):
            result_state[i] = self.state[state.state.shape[0]*i + count]/state.state[count]
        return QuantumState('array', result_state)
    def show(self):
        print(self.state)

# 酉矩阵操作
# 基本门，以及对基本门的一些操作
class QuantumGate():
    # 初始化只支持几种基本门
    def __init__(self, gate='None', parameter=0):
        if gate == 'H':
            self.operation = 1/(2**0.5) * np.array([[1, 1],[1, -1]], dtype='complex64')
        elif gate == 'CNOT':
            self.operation = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], dtype='complex64')
        elif gate == 'Not':
            self.operation = np.array([[0,1],[1,0]], dtype='complex64')
        elif gate == 'Rs':
            self.operation = np.array([[1, 0],[0, e**((2*pi*1j)/2**(parameter))]], dtype='complex64')
        # 自定义门
        elif gate == 'custom':
            if type(parameter) == type([]):
                self.operation = np.array(parameter, dtype='complex64')
            else:
                self.operation = parameter
        elif gate == 'unit':
            self.operation = np.eye(parameter, dtype='complex64')
        elif gate == 'Toffoli':
            self.operation = np.array([[1,0,0,0,0,0,0,0],
                                        [0,1,0,0,0,0,0,0],
                                        [0,0,1,0,0,0,0,0],
                                        [0,0,0,1,0,0,0,0],
                                        [0,0,0,0,1,0,0,0],
                                        [0,0,0,0,0,1,0,0],
                                        [0,0,0,0,0,0,0,1],
                                        [0,0,0,0,0,0,1,0]], dtype='complex64')
            self.operation = self.operation.T
        # 张量积逆顺序，其实就相当于乘一个置换阵
        # 注意，此处并没有用基本门电路来实现
        # 之后会用基本门电路来实现逆序操作
        elif gate == 'reverse':
            self.reverse(parameter)
        elif gate == 'None':
            self.operation = None

    # 实现|a0...an>到|an...a0>的两个qubit交换互换
    # n>=1
    def exchange0(self, n):
        g1 = QuantumGate(gate='Not').conditioned_left(n)
        g2 = QuantumGate(gate='Not').conditioned_right(n)
        g3 = QuantumGate(gate='Not').conditioned_left(n)
        # 由于是右乘，先g1，再g2，再g3，因此从右到左为g1，g2，g3
        return g3@g2@g1
    
    # 对n个qubit进行倒序操作
    def reverse(self, n):
        g1 = self.exchange0(n-1)
        for i in range(n//2-1):
            g1 = self.exchange0(n-1-2*(i+1)).null_operation_left(i+1).null_operation_right(i+1)@g1
        self.operation = g1

    def kron(self, op1):
        return QuantumGate('array', np.kron(self.operation, op1.operation))

    def size(self):
        return self.operation.shape[0]

    # 用一个第n位左加的qubit来控制门
    # 如果只是左边一位来控制，那么为：
    # 对应[[I1, 0], [0, U]]，U为操作门
    # 其中I1的shape和U相同
    # 如果是左边第n位控制，可以分解成向左扩展n-1个null操作，再左边一位控制
    def conditioned_left(self, n=1):
        g1 = self.null_operation_left(n-1)
        length = g1.operation.shape[0]
        I1 = np.eye(length, dtype='complex64')
        zero = np.zeros([length, length], dtype='complex64')
        return QuantumGate('custom', np.hstack((np.vstack((I1, zero)), np.vstack((zero, g1.operation)))))
    
    # 用一个右加的qubit来控制门
    # 若未加之前对应A，则对A，每隔一行插入一行零（第一行为零），然后每隔一列插入一行零（第一列为零），然后将对角线上插入的零变为1即可
    def conditioned_right(self, n=1):
        g1 = self.null_operation_right(n-1)
        temp = np.array([[0, 0],[0, 1]], dtype='complex64')
        temp = np.kron(g1.operation, temp)
        for i in range(temp.shape[0]//2):
            temp[i*2][i*2]=1    
        return QuantumGate('custom', temp)

    # 左加n个不起作用的qubit对应的门
    # 两个情况对应[[A, 0], [0, A]]（就是I卷积gate）n个情况递归
    def null_operation_left(self, n):
        if n == 0:
            return self
        I = np.eye(2**n, dtype='complex64')
        return QuantumGate('custom', np.kron(I, self.operation))
    
    # 右加n个不起作用的qubit对应的门
    # 两个情况对应gate卷积I，n个情况递归
    def null_operation_right(self, n):
        if n == 0:
            return self
        # n个I卷积
        I = np.eye(2**n, dtype='complex64')
        # 然后gate和I卷积
        return QuantumGate('custom', np.kron(self.operation, I))

    def show(self):
        print(self.operation)

    # 运算符重载
    def __matmul__(self, operation1):
        if type(operation1) == type(self):
            return QuantumGate('custom', self.operation @ operation1.operation)
        else:
            return QuantumState('array', self.operation @ operation1.state)
        # else:
        #     print('NotImplemented')
        #     print(operation1)
# 返回QFT操作对应的矩阵
# 参数为array_length，即要变换的序列长度
def QFT(array_length):
    # 所需线路的个数
    length = int(np.log2(array_length))
    # 当前的矩阵操作
    current = QuantumGate('unit', 2**length)
    for i in range(length):
        # 先左补齐
        sub_current = QuantumGate('H').null_operation_left(i)
        # 再右补齐
        sub_current = sub_current.null_operation_right(length-i-1)
        # print(sub_current)
        for q in range(length - i - 1):
            # print(q)
            temp = QuantumGate('Rs', 2+q)
            # 左补齐
            temp = temp.null_operation_left(i)
            # 右补齐，要经历一次右加null，一次右加contioned，一次右加null
            temp = temp.null_operation_right(q).conditioned_right().null_operation_right(length-q-i-2)
            # 注意这里矩阵乘法的顺序
            # print(temp)
            sub_current = temp@sub_current
        current = sub_current@current
    current = ModEx().reverse(length)@current
    return current

# 逆QFT
# 为QFT的mirror circuit
# 这里不再用基本门来实现
# QFT的共轭转置便是逆QFT
def rev_QFT(array_length):
    qft = QFT(array_length)
    qft.operation = np.conjugate(qft.operation).T
    return qft

# 利用可逆电路实现black_box操作
class ModEx():
    def __init__(self):
        pass

    # 实现功能见44页
    def and_gate(self):
        return QuantumGate('Toffoli')

    def not_gate(self):
        return QuantumGate('Not')

    def or_gate(self):
        stage1 = QuantumGate('Not').kron(QuantumGate('Not')).null_operation_right(1)
        stage2 = QuantumGate('Toffoli')
        stage3 = QuantumGate('Not').kron(QuantumGate('Not')).kron(QuantumGate('Not'))
        return stage3@stage2@stage1

    # 线路分开功能
    def Fanout(self):
        stage1 = QuantumGate('Not').null_operation_left(1).null_operation_right(1)
        stage2 = QuantumGate('Toffoli')
        return stage2@stage1
    
    # 实现|a0...an>到|an...a0>的两个qubit交换互换
    # n>=1
    def exchange0(self, n):
        g1 = QuantumGate(gate='Not').conditioned_left(n)
        g2 = QuantumGate(gate='Not').conditioned_right(n)
        g3 = QuantumGate(gate='Not').conditioned_left(n)
        # 由于是右乘，先g1，再g2，再g3，因此从右到左为g1，g2，g3
        return g3@g2@g1

    # 实现|b0^n-1>到|b0 0 b1 0...bn-1>的exchange
    # 首先bn_1和最后的交换，然后bn-2和倒数第三个交换...
    def exchange1(self, n):
        g1 = self.exchange0(n-1).null_operation_left(n-1).null_operation_right(1)
        for i in range(n-2):
            g2 = self.exchange0(n-2-i).null_operation_left(n-2-i).null_operation_right(2*i+3)
            g1 = g2@g1
        return g1
    
    # exchange1的逆
    # 是exchange的mirror
    def rev_exchange1(self, n):
        g1 = self.exchange1(n)
        g1.operation = g1.operation.T
        return g1
    
    # 将|bn-1...b1b0 an-1...a1a0>变为|an-1...a1a0 bn-1..b1b0>
    def exchange2(self, n):
        g1 = self.exchange0(n).null_operation_right(n-1)
        for i in range(n-1):
            g1 = self.exchange0(n).null_operation_right(n-1-i-1).null_operation_left(i+1)@g1
        return g1

    # 对n个qubit进行倒序操作
    def reverse(self, n):
        g1 = self.exchange0(n-1)
        for i in range(n//2-1):
            g1 = self.exchange0(n-1-2*(i+1)).null_operation_left(i+1).null_operation_right(i+1)@g1
        return g1

    # len1是加数的bit长度，即n
    # |q> -> |q+a>
    # x和N内嵌在矩阵里面
    # 2**len1 >= N
    # 返回一个2**len1*2**len1的permutation阵，实现加法
    # len1>=2
    def sum_op(self, a, len1, N):
        # 将a化为nbit表示的形式
        a = a%(2**len1)
        a = bin(a)[2:]
        a = '0'*(len1-len(a))+a
        # 权重小的在前
        a = a[::-1]
        # 加一器
        one_p = (QuantumGate('Not').null_operation_right(1))@(QuantumGate('CNOT'))
        # 构建permutation阵
        # 交换阵，将|bn-1...b1b0>变为|b0b1...bn-1>
        g1 = self.reverse(len1).null_operation_right(len1-1)
        # 交换阵，对2*len1-1个qubit进行操作（有len1-1个ancilla bit）
        g1 = self.exchange1(len1)@g1
        
        if a[0] == '1':
            g1 = (one_p.null_operation_right(2*len1-3))@g1
        for count,i in enumerate(a[1:-1]):
            if i == '1':
                g1 = (one_p.null_operation_left((count+1)*2).null_operation_right(2*len1-1-(count+1)*2-2))@g1
            g1 = (one_p.conditioned_left(1).null_operation_left((count+1)*2-1).null_operation_right(2*len1-1-(count+1)*2-2))@g1
        if a[-1] == '1':
            g1 = QuantumGate('Not').null_operation_left(2*len1-2)@g1
        g1 = QuantumGate('CNOT').null_operation_left(2*len1-3)@g1

        # 权重大的在前
        a = a[::-1]
        # 两种情况
        # 情况1，ai为1
        temp1 = QuantumGate('CNOT').null_operation_right(1)
        temp1 = QuantumGate('Not').null_operation_left(1).null_operation_right(1)@temp1
        temp1 = QuantumGate('CNOT').null_operation_left(1)@temp1
        temp1 = QuantumGate('Not').conditioned_left(2)@temp1
        temp1 = QuantumGate('Toffoli')@temp1
        temp1 = QuantumGate('CNOT').null_operation_right(1)@temp1
        temp1 = QuantumGate('Not').null_operation_left(1).null_operation_right(1)@temp1
        # 情况2，ai为0
        temp2 = QuantumGate('CNOT').null_operation_right(1)
        temp2 = QuantumGate('Toffoli')@temp2
        temp2 = QuantumGate('CNOT').null_operation_right(1)@temp2
        for count,i in enumerate(a[1:-1]):
            if i=='1':
                g1 = temp1.null_operation_right(2*count+1).null_operation_left(2*len1-1-3-2*count-1)@g1
            else:
                g1 = temp2.null_operation_right(2*count+1).null_operation_left(2*len1-1-3-2*count-1)@g1
        if a[-1] == '1':
            temp = QuantumGate('Not').null_operation_right(1)
            temp = QuantumGate('CNOT')@temp
            temp = QuantumGate('Not').null_operation_right(1)@temp
            g1 = temp.null_operation_right(2*len1-3)@g1
        g1 = self.rev_exchange1(len1)@g1
        g1 = self.reverse(len1).null_operation_right(len1-1)@g1
        # 此时g1实现了|b 0^n>到|(b+a) 0^n-1>的映射
        # g1.show()
        return g1

    # 由于时间关系，只构造加法运算，比较和模加运算不再构造
    # 该函数直接返回模加运算的期望矩阵
    def mod_sum_op(self, a, len1, N):
        a = a%N
        g1 = np.zeros([2**len1, 2**len1], dtype='complex64')
        for i in range(N):
            g1[(i+a)%N][i]=1
        for i in range(N, 2**len1):
            g1[i][i]=1
        return QuantumGate('custom', g1)

    def EX_GCD(self, a,b,arr): #扩展欧几里得
        if b == 0:
            arr[0] = 1
            arr[1] = 0
            return a
        g = self.EX_GCD(b, a % b, arr)
        t = arr[0]
        arr[0] = arr[1]
        arr[1] = t - int(a / b) * arr[1]
        return g

    def ModReverse(self, a, n): #ax=1(mod n) 求a模n的乘法逆x
        arr = [0,1,]
        gcd = self.EX_GCD(a,n,arr)
        if gcd == 1:
            return (arr[0] % n + n) % n
        else:
            return -1

    # 乘c操作(mod N)
    # |q> -> |(qc)modN>
    # x和N内嵌在矩阵里面
    def mod_mul_op(self, c, len1, N):
        # 用reverse使权值从小到大排列
        g1 = self.reverse(len1).null_operation_right(len1)
        g1 = self.mod_sum_op(c, len1, N).conditioned_left(len1)@g1
        for i in range(len1-1):
            g1 = self.mod_sum_op(c*2**(i+1), len1, N).conditioned_left(len1-i-1).null_operation_left(i+1)@g1
        # print(g1.operation.shape)
        # print(g1.operation.sum())
        # 加法器的输入和输出，从上到下，从大到小
        # 而控制线路的输入和输出，从上到下，从小到大
        # 因此二者皆需要倒序
        g1 = self.reverse(2*len1)@g1
        # 用扩展欧几里得算法求c乘N的逆
        c_rev = self.ModReverse(c, N)
        for i in range(len1):
            g1 = self.mod_sum_op(-c_rev*2**i, len1, N).conditioned_left(len1-i).null_operation_left(i)@g1
        g1 = self.reverse(len1).null_operation_right(len1)@g1
        # 提取乘法阵，即忽略ancilla bit下的矩阵，方法类比加法
        g1 = g1.operation[::2**(len1),::2**(len1)]
        g1 = QuantumGate('custom', g1)
        # 此时g1输入和输出的权值都是从大到小
        return g1


    # r次方操作(mod N)
    # |r>|1> -> |r>|x^r mod N>
    # x和N内嵌在矩阵里面
    # len1是|r>的长度，len2是|1>的长度
    def mod_exp_op(self, x, len1, len2, N):
        g1 = self.reverse(len1).null_operation_right(len2)
        temp = x
        for i in range(len1):
            g1 = self.mod_mul_op(temp, len2, N).conditioned_left(len1-i).null_operation_left(i)@g1
            temp = (temp*temp)%(N)
        g1 = self.reverse(len1).null_operation_right(len2)@g1
        return g1

# 经典方法寻找周期
# N是odd且not a prime power
# x和N互质
def period_finding(x, N):
    count = 1
    temp = x
    while temp != 1:
        count += 1
        temp = (temp*x)%N
    return count

# Shor方法寻找周期
# 求x^a == x^b mod N(a不等于b)时的(b-a)，(b-a)便是周期
class Shor_period_finding():
    # func为待求周期的函数
    def __init__(self, x, N):
        self.x = x
        self.N = N
        self.l = int(log2(N*N))+1
        self.q = 2**self.l
        self.n = ceil(log2(N))
        print(self.x, self.l, self.n)
        # QFT需要的内存（单位GB），很容易爆内存
        if 128*4**(self.l+self.n)/8589934592 > 1:
            return
    # balck_box，将|a>|0^n>变为|a>|f(a)>，black_box有专门的量子电路实现。
    # 注意，这里的block_box，其实也对应一个置换矩阵
    # |a>|0^n>，用张量积运算后得到的向量形式，应该是[a1,0,0,...,0, a2,0,0,...,0,..., an,0,0,...,0]省略处有2^n-1个零
    # 而|a>|f(a)>，用张量积运算后得到的向量形式，a1会和它之后的2^n-1中的一个零做交换，a2也是如此
    # 这样其实就是对应一个置换矩阵
    # 这里我们先用硬编码的方式实现，之后会再设计量子电路来实现它
    # 它体现了量子计算的并行性，是加速的根本。
    def black_box(self, state):
        exchange = np.zeros([2**(self.l+self.n), 2**(self.l+self.n)])
        for i in range(self.q):
            val = quick_mod(self.x, i, self.N)
            exchange[i*(2**self.n)+val, i*(2**self.n)] = 1
        # 严格来说，这里的exchange并不算是置换阵，因为我们只是对特定的元素进行位置交换，因此该置换阵并不是每行每列都有元素1
        return QuantumGate('custom', exchange)@(state)

    # 利用连分式的原理寻找逼近的分式
    # 返回分式的分子和分母
    def CF_approximate(self, x):
        a0 = int(x)
        p0 = a0
        q0 = 1
        x1 = 1/(x-a0)
        a1 = int(x1)
        p1 = a1*a0+1
        q1 = a1
        p2 = p1
        q2 = q1
        while q2 <= self.N:
            # 此时相等
            if x1 == a1:
                return p1, q1
            x2 = 1/(x1-a1)
            a2 = int(x2)
            p2 = a2*p1 + p0
            q2 = a2*q1 + q0
            x1 = x2
            a1 = a2
            p0 = p1
            p1 = p2
            q0 = q1
            q1 = q2
        # p0/q0就是在保证分母小于等于N的情况下对数x的最佳逼近
        return p0, q0

    def period_finding(self):
        flag = True
        while flag:
            # 以l+n个值为0的qubit的叠加为初始状态
            state = QuantumState('bit', '0'*self.l+'0'*(self.n-1)+'1')
            operation1 = QFT(2**self.l).null_operation_right(self.n)
            # 对前q个qubit进行QFT变换，（or just l Hadamard gates）
            state = operation1@state
            # black_box操作
            operation1 = ModEx().mod_exp_op(self.x, self.l, self.n, self.N)
            state = operation1@state
            # 应用逆QFT
            operation1 = rev_QFT(2**self.l).null_operation_right(self.n)
            state = operation1@state
            # 测量前l个qubit，会坍缩为~s/r（temp1）
            _,temp1 = state.observe([i for i in range(self.l)])
            # 观测坍缩后，看坍缩的b值
            # r为q的因子的情况
                # b/q的值，会等可能取0,1/r，2/r，3/r,...,(r-1)/r
                # 每次取值，得到b/q，当b/q对应的c/r中c和r正好互质时，b/q化成最简形式的分母，就是r值
                # 但是由于会等可能地取0,1/r，2/r，3/r,...,(r-1)/r这些值，而小于r的有r/loglogr个值与r互质
                # 因此有1/loglogr的概率取到互质的值
                # 即尝试loglogr次后，以数学期望1的概率，得到r值
                # 对每次尝试得到的r值，我们还需要检验，即测试x^r mod N是否为1
                # 根据快速模运算即可以检测，这里是用的的经典算法，复杂度是logr
                # 因此可以用logrloglogr的复杂度找到一个因子
                # 由于r<=N，因此复杂度的上界为logNloglogN
                # 也就是说，对于任一个数N，可以用logNloglogN的复杂度找到它的一个因子
                # 更准确地说，复杂度应为MlogNloglogN，M为每次量子操作的复杂度（因为每次都是从0态开始，经过QFT等一系列过程，这也需要复杂度）
            # r不为q的因子
                # 此时，若观测后坍缩，有很高的概率下式成立:|b/q-c/r|<=1/2*q
                # 而在[b/q-1/2q, b/q+1/2q]内，有且只有一个分母小于等于N的真分数，就是c/r  
            count = 0
            while temp1.state[count] == 0:
                count += 1
            # count就是~s/r * 2**l
            # 一直以(m/q)**0.5的概率，b值取零，b值取零对等式求解没有意义
            print(count)
            if count == 0:
                continue
            # 利用连分式的方式进行逼近
            c,r = self.CF_approximate(count/self.q)
            print(c, r)
            # 大概率|b/q-c/r|<=1/2*q成立，该式成立的条件下若c和r互质，那么r就求出了
            if quick_mod(self.x, r, self.N) == 1:
                flag = False
        return r


# 利用period_finding来寻找到N的一个质因子
def period_finding_to_factoring(N, classic):
    if N % 2 == 0:
        return [2]
    # 接下来N不是odd，也不是prime_power
    while 1:
        # 任选x属于{2,...,N-1}
        x = random.randint(2,N-1)
        # 如果不互质，那么最大公因子就是N的一个因子，返回即可
        temp = my_gcd(N, x)
        if temp != 1:
            return [temp]
        # 如果互质，那么x模N乘构成一个群，记x的阶为r，则r有>=0.5的概率为偶数（如果为奇数重新抽取）
        # 经典方法找周期
        if classic:
            period = period_finding(x, N)
        else:
            SPF = Shor_period_finding(x, N)
            period = SPF.period_finding()
        if period % 2 != 0:
            continue
        # a和b不一定为质数，但一定为N的因子
        # print(x,period)
        a = quick_exponential_operation(x, period//2)-1
        b = quick_exponential_operation(x, period//2)+1
        # not multiples of N
        if a%N == 0 or b%N==0:
            continue
        a = my_gcd(a, N)
        b = my_gcd(b, N)
        return [a,b]
        
# 返回N的质因子组成的列表
def Shor_factoring(N, classic = False):
    a,b = is_prime_power(N)
    if a:
        return [b[0] for i in range(b[1])]
    result = []
    # 存储未分解为质数的因子
    my_queue = [N]
    while len(my_queue) != 0:
        temp = my_queue[0]
        my_queue.pop()
        # prime_power的形式可以立即返回
        a,b = is_prime_power(temp)
        if a:
            result+=[b[0] for i in range(b[1])]
            temp //= b[0]**b[1]
        else:
            factor = period_finding_to_factoring(temp, classic)
            # print(factor)
            for i in factor:
                if is_prime(i):
                    result.append(i)
                else:
                    my_queue.append(i)
                temp //= i
        if temp != 1:
            if is_prime(temp):
                result.append(temp)
            else:
                my_queue.append(temp)
    return result


if __name__ == '__main__':
    print(15,'=',Shor_factoring(15))
    # for i in range(2,32):
    #     print(i,'=',Shor_factoring(i))

    # # |r>|1> -> |r>|x^r mod N>
    # # # 模乘运算功能，|r>|1> -> |r>|3^r mod 4>
    # # r有三个qubit来存储
    # a = ModEx().mod_mul_op(3,3,5)
    # print(a.operation)
    # print(a.operation.shape)
    # print(a.operation.sum())
    # # b为2，最后一个bit为ancilla bit
    # for i in range(8):
    #     b = QuantumState('bit', '0'*(3-len(bin(i)[2:]))+bin(i)[2:]+'01')
    #     (a@b).show()

 
