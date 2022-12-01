# KuramotoUKF

### Statistical inference from a coupled-oscillator model
Parameter recovery, latent variable tracking and model comparison of a Kuramoto oscillator network by Unscented Kalman Filter encapsulated in a Markov Chain Monte Carlo algorithm for model comparisons.

The Kuramoto oscillator is a mathematical model commonly used to study synchronisation. As such the model fits well with prominent neuroscience concepts and has potential clinical implications in the study of, for example, pathological Tremor. Kuramoto's original model posed homogenous coupling between oscillators, and while this can be trivially extended to incorporate heterogeneous relationships the analytic solutions quickly break down and a more computational approach is required. The aim of this project was to derive/infer the activity of oscillators as they evolve over time and spontaneously couple/decouple resulting in sporadic bouts of tremor. The generative model is the average activity over all oscillators. We suggest that this representation corresponds to the physical hand movement of a patient with a tremulous condition, such as Parkinson's disease or Essential Tremor, and reflects groups (or pockets) of synchronised oscillators that have been observed in human microelectrode recordings in the Thalamus. Such inference would allow us to experimentally validate the impact of desynchronising pulses on this network (such as delivered through peripheral nerve stimulation, which has previously been trialled for tremor suppression), and derive the optimal conditions to restore homeostasis. This multithreaded code was written to provide high-performance computation for a specific neural model capable of complex interactions.

Implemented techniques:
* Generalised Kuramoto oscillators [[http://en.wikipedia.org/wiki/Kuramoto_model](en.wikipedia.org/wiki/Kuramoto_model)]
* Unscented Kalman Filter (UKF) [[groups.seas.harvard.edu/courses/cs281/papers/unscented.pdf](groups.seas.harvard.edu/courses/cs281/papers/unscented.pdf)]
* Markov Chain Monte Carlo (MCMC) algorithm (Metropolis-Hastings)
* Particle Swarm Optimisation (PSO)
* Gradient Descent (Vanilla, Adam, RMSProp)
* Custom high-performance Maths library

Original repository with commit history: [[http://bitbucket.org/jsbrittain/kuramotoukf/src/master/](bitbucket.org/jsbrittain/kuramotoukf/src/master/)].
