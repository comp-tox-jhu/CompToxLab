
!!! example "Prerequisites"
    - [Deep Learning Setup](./00_setup.md) - Create R project, setup workspace, and download data

## 1. Tensors: Creation and Operations
### Creating Tensors
Tensors are multi-dimensional arrays and the core of PyTorch.

```{python}
scalar = torch.tensor(7)  # Scalar
vector = torch.tensor([1, 2, 3])  # Vector
matrix = torch.tensor([[1, 2], [3, 4]])  # Matrix
tensor = torch.tensor([[[1, 2, 3], [4, 5, 6], [7, 8, 9]]])  # 3D Tensor
```

!!! example "Math Attack!"
    - Scalar: A single number $s \in \mathbb{R}$.
    - Vector: A 1D array $\mathbf{v} \in \mathbb{R}^n$.
    - Matrix: A 2D array $\mathbf{M} \in \mathbb{R}^{m \times n}$.

!!! info "Tips and Tricks"
    - Dimensions as brackets:
    - `[]` → 1D (vector)
    - `[[]]` → 2D (matrix)
    - `[[[]]]` → 3D


## Basic Operations
Perform element-wise operations on tensors.

```{python}
tensor = torch.tensor([1, 2, 3])
print(tensor + 10)  # Addition
print(tensor * 10)  # Multiplication
```

!!! example "Math Attack!"
    - Element-wise addition: $\mathbf{A} + c$, where $c$ is added to each element in $\mathbf{A}$.
    - Element-wise multiplication: $\mathbf{A} \times c$, where $c$ multiplies each element of $\mathbf{A}$.

!!! info "Tips and Tricks"
    - Reassign to modify: Operations don’t change the tensor unless you reassign the result back to some variable.


### Matrix Multiplication
Matrix multiplication combines data, requiring matching inner dimensions.

```{python}
A = torch.tensor([[1, 2], [3, 4]])
B = torch.tensor([[5, 6], [7, 8]])
result = torch.matmul(A, B)
print(result)
```

!!! example "Math Attack!"
    - Matrix multiplication: For $\mathbf{A} \in \mathbb{R}^{m \times n}$ and $\mathbf{B} \in \mathbb{R}^{n \times p}$, the matrix product $\mathbf{C} = \mathbf{A} \cdot \mathbf{B}$ is given by:
    - $c_{ij} = \sum_{k=1}^{n} a_{ik} \cdot b_{kj}$

​
 
!!! info "Tips and Tricks"
    - Shape matching rule: For multiplication, the inner dimensions must match:
    - $A(m \times n) \cdot B(n \times p) \rightarrow C(m \times p)$
    
## 2. Tensor Manipulation: Reshaping, Stacking, and Indexing
### Reshape and Squeeze
Change tensor shapes or remove/add dimensions of size 1.

```{python}
x = torch.arange(1, 10).reshape(3, 3)  # Reshape to 3x3
x_squeezed = x.squeeze()  # Remove dimensions of size 1
x_unsqueezed = x.unsqueeze(0)  # Add a dimension
```

!!! Math Attack!
    - Reshape: Adjusts tensor dimensions while keeping the same number of elements, i.e., $m \times n = p \times q$.
    - Squeeze: Removes dimensions with size 1.
    
!!! info "Tips and Tricks"
    - Squeeze and unsqueeze: Squeeze removes "extra" dimensions, unsqueeze adds them back.

### Stacking Tensors
Stack multiple tensors along a new dimension.

```{python}
x_stacked = torch.stack([x, x], dim=0)
print(x_stacked)
```

!!! info "Tips and Tricks"
    - Stacking is like a sandwich: You’re adding a new layer (dimension) on top of existing tensors.

### Indexing
Select specific elements in tensors.

```{python}
x = torch.tensor([[1, 2, 3], [4, 5, 6]])
print(x[:, 1])  # Second column
print(x[0, 0])  # Element at (0, 0)
```

!!! info "Tips and Tricks"
    - Indexing is slicing: You can slice from any dimension of the tensor, similar to cutting bread slices.

## 3. Tensor Aggregation and Statistics
### Aggregation Functions
Use functions like min, max, mean, and sum to reduce tensors.

```{python}
x = torch.tensor([1, 2, 3, 4])
print(x.min(), x.max(), x.mean(), x.sum())
```

!!! example "Math Attack!"
    - Min: $\text{min}(x)$ returns the smallest element.
    - Max: $\text{max}(x)$ returns the largest element.
    - Mean: $\text{mean}(x) = \frac{1}{n} \sum_{i=1}^{n} x_i$, where $n$ is the number of elements.
    - Sum: $\text{sum}(x) = \sum_{i=1}^{n} x_i$, where $n$ is the number of elements.

### Positional Min/Max
Get the index of the minimum and maximum values.

```{python}
x = torch.tensor([10, 20, 30])
print(x.argmax(), x.argmin())  # Index of max and min
```

## 4. Working with GPUs
### Check GPU Availability
Check if a GPU is available for faster computation.

```{python}
torch.cuda.is_available()  # True if GPU is available
```

### Moving Tensors to GPU
Transfer tensors to the GPU for faster calculations.

```{python}
device = 'cuda' if torch.cuda.is_available() else 'cpu'
tensor = torch.tensor([1, 2, 3]).to(device)
print(tensor)
```

### Moving Back to CPU
Move tensors back to the CPU for further processing.

```{python}
tensor_cpu = tensor.to('cpu')
```
    
## 5. Randomness and Reproducibility
### Set Random Seed
Ensure reproducible results by setting a random seed.

```{python}
torch.manual_seed(42)
random_tensor = torch.rand(3, 4)
print(random_tensor)
```

## 6. Tensor and NumPy Integration
Convert Between NumPy and PyTorch

```{python}
# NumPy to PyTorch
import numpy as np
np_array = np.array([1, 2, 3])
tensor = torch.from_numpy(np_array)
print(tensor)
```

```{python}
# PyTorch to NumPy
numpy_array = tensor.numpy()
print(numpy_array)
```

