#FAQ

##### Q: Why was covar indexing changed?

A: Since version 0.4.4, all indexing is now 1-based instead of 0 based. It turnes out not everyone is a python programmer and most people actually prefer 1-based indeces (go figure. =-) ). The 0-based indexing was causing a lot of confusion, and as such, was changed. Also not that even better than that, you can now provide a file with actual covariate names, which should save you a lot of time further downstream.
