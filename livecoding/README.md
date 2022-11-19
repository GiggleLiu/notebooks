## Generate cast files
```bash
julia -e 'using AsciinemaGenerator; cast_file("einsum.jl"; output_file="einsum.cast", mod=Main, prompt_delay=0.2)'
```