
!!! example "Prerequisites"
    - [Deep Learning Setup](./00_setup.md) : Setup workspace and download python libraries

**Learning Objectives**

1. [What is a Tensor?](#what-is-a-tensor)
2. [Basic Operations](#basic-operations)
3. [Matrix Multiplication](#matrix-multiplication)
4. [Tensor Manipulation](#tensor-manipulation)
5. [Accessing Elements in a Tensor](#accessing-elements-in-a-tensor)
6. [Tensor Statistics](#tensor-statistics)
7. [Working with GPUs](#working-with-gpus)
8. [Random Values and Reproducibility](#random-values-and-reproducibility)
9. [Tensor and NumPy Integration](#tensor-and-numpy-integration)

## What is a Tensor?

Machine learning is all about manipulating numbers and to do that we have different ways of storing those numbers, typically in structures called tensors. Tensors can be a single number (scalar), a list or vector of numbers (1D tensor), a matrix (2D tensor), a list of matrices (3D tensor), so on and so forth:

!!! info "What is a Tensor?"
    <figure markdown="span">
      ![](img/tensors.png){ width="500" }
    </figure>

Now let's see how we make a tensor in PyTorch!

```{python}
scalar = torch.tensor(7)  # Scalar
vector = torch.tensor([1, 2, 3])  # Vector
matrix = torch.tensor([[1, 2], [3, 4]])  # Matrix
tensor = torch.tensor([[[1, 2, 3], [4, 5, 6], [7, 8, 9]]])  # 3D Tensor
```

Machine learning papers always include a lot of math jargon, but don't be afraid! Let's go through a few of the symbols for tensors! 

- **Scalar**: A single number $s \in \mathbb{R}$.
- **Vector**: A 1D array or vector is shown as: $\mathbf{v} \in \mathbb{R}^n$.
- **Matrix**: A 2D array or matrix with m rows and n columns is shown as: $\mathbf{M} \in \mathbb{R}^{m \times n}$.
- **3D Tensor**: a list of o matrices with m rows and n columns is shown as: $\mathbf{M} \in \mathbb{R}^{o \times m \times n}$.

Now in python these tensors are shown using different numbers of brackets:

!!! info "Tips and Tricks"
    - Dimensions as brackets:
    - `[]` → 1D (vector)
    - `[[]]` → 2D (matrix)
    - `[[[]]]` → 3D


## Basic Operations

Again, tensors are just numbers and we can perform basic operations on them, like addition, subtraction, multiplication and division. Let's go through how we would do that!

```{python}
tensor = torch.tensor([1, 2, 3])
print(tensor + 10)  # Addition
print(tensor * 10)  # Multiplication
```

!!! info "output"
    ```{python}
    tensor([11, 12, 13])
    tensor([10, 20, 30])
    ```

Now in math speak we represent the operations above (and others like division and subtraction) like so: 

- **Element-wise addition**: $\mathbf{A} + c$, where $c$ is added to each element in $\mathbf{A}$.
- **Element-wise multiplication**: $\mathbf{A} \times c$, where $c$ multiplies each element of $\mathbf{A}$.

!!! info "Tips and Tricks"
    - Reassign to modify: Operations don’t change the tensor unless you reassign the result back to some variable.


## Matrix Multiplication

A lot of machine learning relies heavily on matrix operations and often times multiplying matrices together. When we multiple matrices we need to match the inner dimensions and then we end up with a matrix with the outer dimensions:

!!! success "This Works!"
    $A(4 \times 2) \cdot B(2 \times 3) \rightarrow C(4 \times 3)$
    This works because the inner dimensions are 2 and 2

!!! failure "This Doesn't Work!"
    $A(4 \times 3) \cdot B(2 \times 3) \rightarrow X$
    This does not work because the inner dimensions are 3 and 2. They are different!

Let's do this in PyTorch!

```{python}
# Define A as a 2x3 tensor
A = torch.tensor([[1, 2, 3], [4, 5, 6]])
# Define B as a 3x2 tensor
B = torch.tensor([[7, 8], [9, 10], [11, 12]])
# Perform matrix multiplication
result = torch.matmul(A, B)
print(result)
```

!!! info "output"
    ```{python}
    tensor([[ 58,  64],
            [139, 154]])
    ```

!!! tip "Tips and Tricks"
    - Instead of writing out `torch.matmul(A,B)` you can use the shorthand `A @ B` for the same result!

Now what would this look like in math speak?

- Matrix multiplication:
- $\mathbf{A} \in \mathbb{R}^{m \times n}$
- $\mathbf{B} \in \mathbb{R}^{n \times p}$
- the matrix product is:
- $\mathbf{C} = \mathbf{A} \cdot \mathbf{B}$
- and each value is calculated by
- $c_{ij} = \sum_{k=1}^{n} a_{ik} \cdot b_{kj}$

Now that's a lot, let's just look at a visualization of how you multiply matrices:
    
<video controls>
<source src="../img/mat_mult.mp4" type="video/mp4">
</video>

## Tensor Manipulation

Ok so when working with tensors and especially performing matrix operations, we will often be trying to match dimensions. We can do that using a number of different functions in PyTorch. Let's start out by making a vector:

```{py}
x = torch.tensor([1,2,3,4,5,6])
x
```

!!! info "output"
    ```
    tensor([1, 2, 3, 4, 5, 6])
    ```

Great, now let's add an extra dimension so that it no longer a vector but one row of a matrix. We can do this by using the `view` or `reshape`function:

```{py}
x.view(1,6),x.reshape(1,6)
```

!!! info "output"
    ```
    (tensor([[1, 2, 3, 4, 5, 6]]), tensor([[1, 2, 3, 4, 5, 6]]))
    ```
Now what about making it a column?

```{py}
x.view(6,1),x.reshape(6,1)
```

!!! info "output"
    ```
    (tensor([[1],
             [2],
             [3],
             [4],
             [5],
             [6]]),
     tensor([[1],
             [2],
             [3],
             [4],
             [5],
             [6]]))
    ```
Now what if we wanted to convert it back to just that list of numbers?

```{py}
x.view(1,6).squeeze()
```

!!! info "output"
    ```
    tensor([1, 2, 3, 4, 5, 6])
    ```
Hmmm, I want to have my column back, let's reverse with unsqueeze! Here we set dim=1 so that it unsqueezes back to a column, dim=0 would make it a row:

```{py}
x.view(1,6).squeeze().unsqueeze(dim=1)
```

!!! info "output"
    ```
    (tensor([[1],
             [2],
             [3],
             [4],
             [5],
             [6]])
    ```
Now what if we needed to stack our tensors? Well the `stack` function will do just that. Here we stack by row:

```{py}
torch.stack([x,x],dim=0)
```

!!! info "output"
    ```
    tensor([[1, 2, 3, 4, 5, 6],
          [1, 2, 3, 4, 5, 6]])
    ```
Let's stack by column as well!   

```{py}
torch.stack([x,x],dim=1)
```

!!! info "output"
    ```
    tensor([[1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
            [5, 5],
            [6, 6]])
    ```
## Accessing Elements in a Tensor

To get specfic elements in a tensor we should start by taking an example tensor:

```{py}
one_mat = torch.tensor([[[1,2,3],
                         [4,5,6],
                         [7,8,9]]])
one_mat
```

!!! info "output"
    ```
    tensor([[[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]]])
    ```

Let's see it's shape:

```{py}
one_mat.shape
```

!!! info "output"
    ```
    torch.Size([1, 3, 3])
    ```

Now let's see what is the first element in the tensor:

```{py}
one_mat[0]
```

!!! info "output"
    ```
    tensor([[1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]])
    ```
Hey look at that the first element is the matrix, let's look at the first element of this first element:

```{py}
one_mat[0][0]
```

!!! info "output"
    ```
    tensor([1, 2, 3])
    ```

Now how about the first element of this element:

```{py}
one_mat[0][0][0]
```

!!! info "output"
    ```
    tensor(1)
    ```
    
## Tensor Statistics

Now often times we may want to summarize elements in our tensors. What is the maximum or minimum value? How about the sum? The mean? We can do this using the following functions:

```{python}
# here we specify that the data type is a float so that we can get these summary statistics!
x = torch.tensor([1, 2, 3, 4], dtype=torch.float)
print(x.min(), x.max(), x.mean(), x.sum())
```

!!! info "output"
    ```
    tensor(1.) tensor(4.) tensor(2.5000) tensor(10.)
    ```

- **Minimum**: $\text{min}(x)$ returns the smallest element.
- **Maximum**: $\text{max}(x)$ returns the largest element.
- **Mean**: $\text{mean}(x) = \frac{1}{n} \sum_{i=1}^{n} x_i$, where $n$ is the number of elements.
- **Sum**: $\text{sum}(x) = \sum_{i=1}^{n} x_i$, where $n$ is the number of elements.


It may be important to grab _where_ in the tensor our minimum or maximum is as well:

```{python}
x = torch.tensor([10, 20, 30])
print(x.argmax(), x.argmin())  # Index of max and min
```

!!! info "output"
    ```
    tensor(2) tensor(0)
    ```

## Working with GPUs

Machine learning can be hard and often needs computer power to make it happen. One way it can make that happen is with graphic processing units or GPUs. These are hardware that are great at handling certain tasks like matrix math really fast. So it can be useful to move our tensors to a GPU so that our AI models can run faster! But first let's check to see if they are available:

```{python}
torch.cuda.is_available()  # True if GPU is available
```

!!! info "output"
    ```
    True
    ```

Now to move a tensor to GPU:

```{python}
device = 'cuda' if torch.cuda.is_available() else 'cpu'
tensor = torch.tensor([1, 2, 3]).to(device)
print(tensor)
```
To move it back to the central processing unit or CPU, just do the following:

```{python}
tensor_cpu = tensor.to('cpu')
```
    
## Random Values and Reproducibility

Now we often times initialize model values with random numbers to ensure those numbers are consistent, you can set a "seed" number to ensure consistency:

```{python}
torch.manual_seed(42)
random_tensor = torch.rand(3, 4)
print(random_tensor)
```

!!! info "output"
    ```
    tensor([[0.8823, 0.9150, 0.3829, 0.9593],
            [0.3904, 0.6009, 0.2566, 0.7936],
            [0.9408, 0.1332, 0.9346, 0.5936]])
    ```
    
## Tensor and NumPy Integration

Tensors are special objects designed for complex operations. However, we may want to convert them back to numpy objects so that we can do things like plot model results. First we will show you how to take a numpy array and make a tensor:

```{python}
# NumPy to PyTorch
import numpy as np
np_array = np.array([1, 2, 3])
tensor = torch.from_numpy(np_array)
print(tensor)
```

!!! info "output"
    ```
    tensor([1, 2, 3], dtype=torch.int32)
    ```

Now let's convert that back to numpy!

```{python}
# PyTorch to NumPy
numpy_array = tensor.numpy()
print(numpy_array)
```

!!! info "output"
    ```
    [1, 2, 3]
    ```

**Key Points**

- Tensors store data across different dimensions: scalars (0D), vectors (1D), matrices (2D), and higher.
- Basic element-wise operations can be used to modify tensors (i.e. addition, subtraction, multiplication, division, etc.).
- For matrix multiplication, the inner dimensions must match!
- Tensor shapes can be manipulated which can be necessary for things like matrix multiplication.
- Tensor elements are accessed by the bracket they are in
- GPUs can be used for faster computation.
- Set random seeds with torch.manual_seed() for reproducibility.
- Converting between PyTorch tensors and NumPy arrays can be useful when plotting model metrics
