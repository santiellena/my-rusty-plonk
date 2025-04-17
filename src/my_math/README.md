## About this folder...

I wrote this files to learn from zero the math of zkSNARKs, specifically the Plonk zkSNARK under the KZG Polynomial Commitment Scheme. Note that the operations are done in 64 bits numbers which is not useful for the actual prime numbers and elliptic curve orders that are used in real schemes (typically BN254).

I found it easier to start like this and then move on into bigger fields because everything was more intuitive and I'm kind of new to Rust, so I didn't want to go and work with integrations right away.

Hope this set of files are useful and help someone in his/her learning path.

Some resources that could be helpful to get this whole thing:

- [ZK Hack MoonMath Study Group](https://zkhack.dev/zk-study-group-moonmath-manual/)
- [ZK Book from RareSkills](https://www.rareskills.io/zk-book)
- Electisec Fellowship Recordings from [1st Cohort](https://www.youtube.com/playlist?list=PLeUIc0UZxuuF8_ueHNt1TuEyNhcsmzu_g) and [2nd Cohort](https://www.youtube.com/playlist?list=PLeUIc0UZxuuEfY9GbIrYaa9tKXT75d-MB)
- [The Discrete Logarithm Problem by Khan Academy](https://youtu.be/SL7J8hPKEWY?si=01cU525dzrjXvU9r)
- [Euclidean Division by LambdaClass](https://www.youtube.com/watch?v=5W4ghK7dWuI&t=1s)
- [Lagrange Interpolation by LambdaClass](https://www.youtube.com/watch?v=REnFOKo9gXs)
- [Exploring Elliptic Curve Pairings by Vitalik Buterin](https://medium.com/@VitalikButerin/exploring-elliptic-curve-pairings-c73c1864e627)
- [Ethereum Foundation bn128 Field Elements Implementation](https://github.com/ethereum/py_pairing/blob/master/py_ecc/bn128/bn128_field_elements.py)
- [Plonk By Hand Series](https://research.metastate.dev/plonk-by-hand-part-1/). Special love for this one ❤.

For the real and useful Plonk, see [this folder](../math_with_ark/). There you will find the Plonk specific resources! However, getting the math from here, is a pre-requisite, so do your homework first!