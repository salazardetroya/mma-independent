# project.py
import flow
from flow import FlowProject
import itertools
import platform
import sys

sys.settrace

plt = platform.system()


def label_text(nref, inner_product):
    return f"{nref}_{inner_product}"


@FlowProject.label
def check_finished(job):
    check = (
        job.isfile(f"finished_{label_text(job.sp.nref, job.sp.inner_product)}.txt")
        and job.isfile(f"output_{label_text(job.sp.nref, job.sp.inner_product)}.txt")
    )
    return check


@FlowProject.operation
@flow.cmd
@FlowProject.post(check_finished)
def launch_opti(job):
    if plt == "Linux":
        simulation = "source /g/g92/miguel/workspace/firedrake_setup.sh \n"
    else:
        simulation = ""

    n_ref_to_n_procs = {1 : 1, 2: 2, 3: 36, 4: 200}

    output = f"{job.ws}/output_{label_text(job.sp.nref, job.sp.inner_product)}.txt"
    if plt == "Linux":
        simulation += f"srun -n {n_ref_to_n_procs[job.sp.nref]} --output={output} "
    simulation += f"python3 compliance_3D.py \
            --output_dir {job.ws} \
            --inner_product {job.sp.inner_product} \
            --nref {job.sp.nref} "
    if plt == "Linux":
        simulation += "\n"
    else:
        simulation += (
            f" &> {job.ws}/output_{label_text(job.sp.nref, job.sp.inner_product)}.txt\n"
        )
    return simulation


if __name__ == "__main__":
    FlowProject().main()
