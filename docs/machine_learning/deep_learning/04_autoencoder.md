
!!! example "Prerequisites"
    - [Deep Learning Setup](./00_setup.md) : Setup workspace and download python libraries

**Learning Objectives**
1. [XXX](#)
2. [XXX](#)

## Introduction

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


So, why bother - how does this contribute to prevent overfitting? Well, by encoding and decoding the input data we tend to "denoise" our data, allowing the model to learn general patterns without memorizing our data. To assess an autoencoder we typically use the mean squared error (MSE) for our loss:

$$
L_{reconstruction} = \frac{1}{N} \sum_{i=1}^N (x_i - \hat{x}_i)
$$

- $L_{reconstruction}$: loss after encoding and decoding
- $N$: number of observations
- $x_i$: true value
- $\hat{x}_i$: predicted value after endcoding and decoding

With a lower $L_{reconstruction}$, the model is doing a better job of reconstructing the original data!


