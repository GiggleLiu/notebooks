## Generate cast files
```bash
julia -e 'using AsciinemaGenerator, InteractiveUtils; cast_file("einsum.jl"; output_file="einsum.cast", mod=Main, tada=true)'
```