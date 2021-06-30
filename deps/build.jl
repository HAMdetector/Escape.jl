using Conda

Conda.add("tensorflow")
Conda.pip_interop(true)
Conda.pip("install", "mhcflurry")

run(`$(Conda.BINDIR)/mhcflurry-downloads fetch`)

run(`$(Conda.BINDIR)/mhcflurry-downloads fetch models_class1_presentation`)