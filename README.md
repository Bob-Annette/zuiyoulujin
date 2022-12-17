Hello，大家好！读研前写过一篇遗传算法的代码，比较简单，算是个入门，当时就有想用它来解决最优路径的问题，上算法导论课时碰巧有听到同学有分享过，但由于自己研究的方向不是这块，就没有再弄，结果今年的华为杯数学建模竞赛F题居然有所涉及，真是用时方恨晚，最近即将毕业，也稍微空闲些了，就再用遗传算法慢慢捡回我的公众号。今天分享的是如何用遗传算法进行最优路径求解问题！这是两年前写的遗传算法，做了一个简单的介绍，感兴趣的小伙伴可以翻看。

# 问题

现在有$n$个地址的坐标，以第一个为起点，途径所有地址，再回到起点，所有地方仅去一次，规划最短路径。

# 思路与Python实现

## 编码

首先先解决编码问题，与上篇文章不同，这次解决的是路径规划问题，有的是一个一个的坐标点，因此我们采用“符号编码”代表这些坐标点，染色体上的编码顺序代表路径顺序。

我们随机生成十组坐标，用作本文的示范：

```python
a = list(range(-5,5))
b = list(range(-4,6))
random.shuffle(a)
random.shuffle(b)
local = list(zip(a,b))
```

| 序号 | 横坐标 | 纵坐标 |
| :--: | :----: | :----: |
|  0   |   -5   |   -2   |
|  1   |   -3   |   4    |
|  2   |   4    |   5    |
|  3   |   -4   |   -3   |
|  4   |   -1   |   1    |
|  5   |   0    |   -1   |
|  6   |   1    |   3    |
|  7   |   2    |   -4   |
|  8   |   -2   |   0    |
|  9   |   3    |   2    |

这些序号就可以用作我们符号编码，例如[0,1,2,3,4,5,6,7,8,9]

## 个体评价

个体评价也就是我们的目标函数，用来区分群体中个体好坏的标准。我们路径规划问题是寻找最短行驶路径，这里用两点之间距离公式进行度量，然后按照染色体上的编码顺序依次累加两点之间的距离：

```python
def dis(start, end):
    # 两点之间距离
    return np.sqrt((end[0]-start[0])**2+(end[1]-start[1])**2)

def fuc(x):
    dis_sum = dis((0,0),local[x[0]])
    for i in range(1,len(x)):
        dis_sum += dis(local[x[i-1]],local[x[i]])
    dis_sum += dis(local[x[-1]],(0,0))
    return dis_sum

def get_fitness(pops):
    return list(map(fuc,pops))
```

*注：本文是以原点作为线路的起点和终点，形成闭环，不同情况要不同设计*

## 选择

选择算子的作用是对个体进行优胜劣汰：从父代群体中选取一些适应度高个体遗传到下一代群体中，这次采用**锦标赛选择策略**。

**锦标赛选择策略**：从种群中随机采样$s$个个体（有放回抽样）进行PK，其中适应度值最优的个体胜出，成为下一代的父代基因，进行$k$轮，得到$k$个优质父代。

这样最差的个体永远不会存活，并且计算简单，不容易陷入局部最优，可以达到更好的求解效果。

```python
def select(pops):
    k = round(np.sqrt(len(pops)+0.25)+0.5)
    fitness = get_fitness(pops)
    father_pops = []
    for i in range(k):
        min_index = np.array(fitness).argmin()
        father_pops.append(pops.pop(min_index))
        fitness.pop(min_index)
    return father_pops
```

## 交叉

接下来就到了整个算法的重头，不同于上一篇文章，单纯地使用两点交叉即可，路径规划中，所有的地点要秉持“不遗漏，不重复”的原则，如果单纯地交叉，会导致地点重复或遗漏。因此就在两点交叉的过程中加一个映射过程，如下图所示：

![两点交叉策略示例](.\figures\两点交叉策略示例.png)

```python
def yinshe(dic, x): # 映射
    while x in dic.keys():
        x = dic[x]
    return x

def jiaocha(Lis1, Lis2): # 两点交叉
    lis1 = Lis1.copy()
    lis2 = Lis2.copy()
    n = len(lis1)
    cross_points_1 = np.random.randint(low=0, high=n-1)
    cross_points_2 = np.random.randint(low=cross_points_1+1, high=n)
    
    yinshe2 = dict(zip(lis1[cross_points_1: cross_points_2],lis2[cross_points_1: cross_points_2]))
    yinshe1 = dict(zip(lis2[cross_points_1: cross_points_2],lis1[cross_points_1: cross_points_2]))

    lis1[cross_points_1: cross_points_2], lis2[cross_points_1: cross_points_2] = lis2[cross_points_1: cross_points_2], lis1[cross_points_1: cross_points_2]
    
    lis1[:cross_points_1] = list(map(lambda x: yinshe(yinshe1, x), lis1[:cross_points_1]))
    lis1[cross_points_2:] = list(map(lambda x: yinshe(yinshe1, x), lis1[cross_points_2:]))
    lis2[:cross_points_1] = list(map(lambda x: yinshe(yinshe2, x), lis2[:cross_points_1]))
    lis2[cross_points_2:] = list(map(lambda x: yinshe(yinshe2, x), lis2[cross_points_2:]))
    return lis1, lis2
```

细心的小伙伴肯定就发现了，为什么映射过程的代码里有一个while循环，这个地方是我在实验的过程中发现的一个细节，如果单靠一次映射，并不能保证所重复的点都映射完，就像上图中的例子，明明$9$对应的是$4$，但是$4$本身也在交叉的基因之中，并没能把9映射到外面，故需要三次映射，即$9 \rightarrow 4 \rightarrow 3 \rightarrow 6$，因此才把外面的$9$替换为了$6$。

## 变异

是对群体中的个体的某些基因座上的基因值作变动，模拟生物在繁殖过程，新产生的染色体中的基因会以一定的概率发生突变。这样的设计可以很好地避免局部最优的情况。

```python
def mutation(pop, MUTATION_RATE=0.003):
    if np.random.rand() < MUTATION_RATE:
        n = len(pop)
        cross_points_1 = np.random.randint(low=0, high=n-1)
        cross_points_2 = np.random.randint(low=cross_points_1+1, high=n)
        pop[cross_points_1], pop[cross_points_2] = pop[cross_points_2], pop[cross_points_1]
    return pop
```

## 整合

有了“个体评价”、“选择”、“交叉”、“变异”这些模块，就可以实现一代代的遗传进化：

```python
def evolution(pops):
    father_pops = select(pops)
    k = len(father_pops)
    new_pops = []
    for i in range(k-1):
        for j in range(i+1,k):
            son1, son2 = jiaocha(father_pops[i], father_pops[j])
            son1 = mutation(son1, MUTATION_RATE=0.3)
            son2 = mutation(son2, MUTATION_RATE=0.3)
            new_pops.append(son1)
            new_pops.append(son2)
    return new_pops
```

## 效果演示

最后只需要把前面的内容整合在一起即可，因为问题比较简单，所以遗传的代数设置50就足够了。

```python
local = [(-5, -2), (-3, 4), (4, 5), (-4, -3), (-1, 1), (0, -1), (1, 3), (2, -4), (-2, 0), (3, 2)]
pops = initpop(N, list(range(10)))
for _ in range(50):
    pops = evolution(pops)
print(min(get_fitness(pops)))
```

然后我们把整个遗传过程可视化出来，效果如所想的一样，最短的路径就是围着转一圈，整个过程感觉还是非常神奇的，如果有想要这个可视化代码的小伙伴，可在文末获得，文章中就不再赘述了。

![演示](.\figures\演示.gif)

# 获得代码

以下是我的个人公众号，本文完整代码已上传，关注公众号回复“**遗传算法最优路径**”，即可获得，谢谢大家支持。

![0](.\figures\0.jpg)
