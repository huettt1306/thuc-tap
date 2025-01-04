import subprocess

def create_slurm_script(script_name, job_name, commands, output_dir):
    """
    Tạo script SLURM từ các lệnh được cung cấp.

    Parameters:
        script_name (str): Đường dẫn để lưu script SLURM.
        job_name (str): Tên job trong SLURM.
        commands (str): Các lệnh thực thi trong SLURM script.
        output_dir (str): Thư mục để lưu log output và error.
    """
    slurm_script = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output={output_dir}/{job_name}.log
#SBATCH --error={output_dir}/{job_name}.err

set -e  # Dừng nếu gặp lỗi

{commands}
"""
    with open(script_name, "w") as file:
        file.write(slurm_script)

def submit_to_slurm(script_path):
    """
    Submit script SLURM lên hệ thống SLURM.

    Parameters:
        script_path (str): Đường dẫn tới file SLURM script cần submit.

    Returns:
        str: Job ID được gán bởi SLURM.
    """
    try:
        result = subprocess.run(["sbatch", script_path], capture_output=True, text=True, check=True)
        output = result.stdout.strip()
        if output.startswith("Submitted batch job"):
            job_id = output.split()[-1]
            return job_id
        else:
            raise ValueError(f"Unexpected output from sbatch: {output}")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job: {e.stderr}")
        raise

def get_job_status(job_id):
    """
    Kiểm tra trạng thái của một job trong SLURM.

    Parameters:
        job_id (str): ID của job cần kiểm tra.

    Returns:
        str: Trạng thái của job (Pending, Running, Completed, Failed, etc.).
    """
    try:
        result = subprocess.run([
            "squeue", "-j", job_id, "-h", "-o", "%T"
        ], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return "UNKNOWN"
    except Exception as e:
        return f"Error: {e}"

def is_job_completed(job_id):
    """
    Kiểm tra xem một job đã hoàn thành hay chưa.

    Parameters:
        job_id (str): ID của job cần kiểm tra.

    Returns:
        bool: True nếu job đã hoàn thành (COMPLETED), False nếu chưa.
    """
    status = get_job_status(job_id)
    return status == "COMPLETED"

def wait_for_job_completion(job_id, interval=60):
    """
    Chờ đến khi job hoàn thành.

    Parameters:
        job_id (str): ID của job cần kiểm tra.
        interval (int): Thời gian (giây) giữa các lần kiểm tra.
    """
    import time

    print(f"Waiting for job {job_id} to complete...")
    while True:
        status = get_job_status(job_id)
        print(f"Job {job_id} status: {status}")
        if status == "COMPLETED":
            print(f"Job {job_id} has completed.")
            break
        elif status in ["FAILED", "CANCELLED"]:
            print(f"Job {job_id} has failed or was cancelled.")
            break
        time.sleep(interval)
