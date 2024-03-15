using PkgTemplates

t = Template(;
    user="mimikhas",
    host="gitlab.cern.ch",
    license="MIT",
    authors="Misha Mikhasenko",
    julia_version=v"1.2",
    dir=".",
    ssh=true,
    plugins=[
        TravisCI(),
        Codecov(),
        AppVeyor(),
    ]
)

generate("JpsiJpsi", t)
isfile(joinpath("C:", "Users", "mikha", "Documents", "JpsiJpsi"))
