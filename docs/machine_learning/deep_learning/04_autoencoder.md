
!!! example "Prerequisites"
    - [Deep Learning Setup](./00_setup.md) : Setup workspace and download python libraries

**Learning Objectives**

1. [XXX](#)
2. [XXX](#)

## Introduction to Autoencoders

Neural networks can be great at learning patterns in data. But the trade off is that the model can be too good, meaning it essentially memorizes all the training data and is not generalizable to other data - in other words the model is overfit:

!!! info "Model Overfitting"
    <figure markdown="span">
      ![](img/overfitting.png){ width="600" }
      <figcaption>Image Credit: [H2O AI](https://h2o.ai/wiki/overfitting/)</figcaption>
    </figure>

So enter autoencoders! Autoencoders take input from a higher dimensional space and _encode_ it in a lower dimensional space, then _decode_ the output of the latent space and reconstructs the data:

!!! info "Autoencoders"
    <figure markdown="span">
      ![](img/autoencoder.png){ width="400" }
      <figcaption></figcaption>
    </figure>


So, why bother - how does this contribute to prevent overfitting? Well, by encoding and decoding the input data we tend to "denoise" our data, allowing the model to learn general patterns without memorizing our data. To assess an autoencoder we use the loss function that makes most sense for our problem. In this tutorial we are going to be predicting smoking status again (a binary variable) so we need the binary cross entropy loss function!

$$
L_{BCE} = - \frac{1}{N} \sum_{i=1}^N [y_i log({\hat{y}_i}) + (1-y_i)log(1- \hat{y}_i)]
$$

- $L_{BCE}$: loss after encoding and decoding
- $N$: number of observations
- $y_i$: true value
- $\hat{y}_i$: predicted value after endcoding and decoding

With a lower $L_{BCE}$, the model is doing a better job of reconstructing the original data! 

## Building an Autoencoder

Let's make ourselves an autoencoder then! We will pull in our glioblastoma data from the classifier example:

```py
# pandas and numpy for data manipulation and plotly for plotting
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# sklearn for model metrics, normalization and data splitting
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

# torch for tensor operations and neural network functions
import torch
from torch.utils.data import DataLoader, TensorDataset
import torch.nn as nn
import torch.optim as optim

# read in data and filter out na values
gbm = pd.read_csv("../data/gbm_data.csv",on_bad_lines='skip')
gbm = gbm[gbm['SMOKING_HISTORY'].notna()]

# make a new column for smoking status
gbm['smoking_status'] = np.where(gbm['SMOKING_HISTORY'].str.contains('non-smoker'), 'non-smoker', 'smoker')

# filter our data for features of interest
gbm_filt = gbm.loc[:,['LINC00159','EFTUD1P1','C20orf202','KRT17P8','RPL7L1P9','smoking_status']]
gbm_filt['smoking_status'] = gbm_filt['smoking_status'].astype('category').cat.codes
```

Now we will split our data into training and test datasets and normalize:

```py
# Split the data into features and outcome variable
X = gbm_filt.drop('smoking_status', axis=1).values
y = gbm_filt['smoking_status'].values

# Normalize the features
scaler = StandardScaler()
X = scaler.fit_transform(X)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=81)

# Convert to PyTorch tensors
X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train, dtype=torch.long)
X_test_tensor = torch.tensor(X_test, dtype=torch.float32)
y_test_tensor = torch.tensor(y_test, dtype=torch.long)
```

Great, let's make that autoencoder now!

```py
class AutoEncoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super(AutoEncoder, self).__init__()
        # encoder 
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, latent_dim)
        )
        # decoder 
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()

        )
    # define the flow through the nodes
    def forward(self, x):
        latent = self.encoder(x)
        reconstructed = self.decoder(latent)
        return reconstructed
```

Here we:

- Create an autoencoder where the encoder takes our initial inputs
- Sends to two hidden layers with 64 nodes and a ReLU activation function
- Squeezes the input to some latent dimension
- Next in the decoder, we take the input from the latent space and send it to two more layers with the original dimensions of the hidden layers - 64
- Then it outputs probabilities for binary classification!
- Finally, we define how is data going to move through our model - it is encoded by our encoder and then decoded by our decoder.

Let's train our autoencoder and see how well it does!

```py
# Correct input dimension
input_dim = X_train_tensor.shape[1]

model = AutoEncoder(input_dim=input_dim, latent_dim=2)
criterion = nn.BCELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

# Convert target tensors to float32
y_train_tensor = y_train_tensor.float()
y_test_tensor = y_test_tensor.float()

loss_vals = []
for epoch in range(1000):
    model.train()

    output = model(X_train_tensor)
    loss = criterion(output.squeeze(), y_train_tensor)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    # Loss for each epoch
    loss_vals.append(loss.item())
```


!!! info "Diagram of Model Workflow"

    ``` mermaid
    flowchart LR
    subgraph Initialization
        A[Initialize Model, Criterion, and Optimizer]
        B[Convert Target Tensors to float32]
    end

    subgraph Training Loop
        C[Loop: 1000 Epochs]
        C -->|Each Epoch| D[Generate Output]
        D --> E[Calculate Loss]
        E --> F[Backpropagate and Update Model]
        F --> G[Record Loss]
        G --> C[Repeat until 1000 Epochs Complete]
    end

    subgraph Output
        H[Output Loss Values]
    end

    A --> B --> C
    G --> H
    ```

You'll see our familiar flow, where we:

- call our model with the our number of input features, and the number of nodes in the latent space
- define our loss function and optimizer
- convert our outcome variables to float values
- loop through our epochs in training mode
- plug and chug our model
- get our loss
- clear our gradients
- backpropogate to update our model weights
- perform gradient descent to find the optimal weights

Now let's see how our loss changes as the number of epochs increases!

```py
# Create the plot
fig = go.Figure()

fig.add_trace(go.Scatter(x=list(range(1000)), y=loss_vals, mode='lines+markers', name='Loss'))

# Add titles and labels
fig.update_layout(
    title='Loss vs. Epoch',
    xaxis_title='Epoch',
    yaxis_title='Loss',
    template='plotly_white')

# Show the plot
fig.show()
```

!!! info "output"
    <figure markdown="span">
      ![](img/autoencoder_loss.png){ width="400" }
      <figcaption></figcaption>
    </figure>

Great! we see that as the number of epochs increases the loss decreases and that this loss seems to level out after epoch 300.

## Variational Autoencoders

A recent improvement on the autoencoder is the variational autoencoder (VAE). So a normal autoencoder is returning a reconstructed value in an effort to avoid overfitting the model. A VAE will turn the latent space into a distribution instead of a vector by returning two things from the latent space:

- **Mean ($\mu$)**: the central point of the distribution.
- **Log variance $log(\sigma^2)$**: describes how spread out or “wide” the distribution should be around the mean.

Now to get data points they essentially sample from a normal distribution defined by the mean and variance output from the VAE and output another data point $z$: 

$z = \mu + \sigma \odot \epsilon$

Where:

- $\mu$ is the mean
- $\sigma$ is the standard variation (which you get by converting the log variance with $exp(0.5⋅log(\sigma^2))$ )
- $\epsilon$ is the sample from $\mathcal{N}(0,1)$ (a normal distribution with a mean of 0 and a standard deviation of 1)
- $\odot$ is the element wise product

What is cool about this is that we can use variational autoencoders to model the underlying distribution. However, the loss may not square up with our 