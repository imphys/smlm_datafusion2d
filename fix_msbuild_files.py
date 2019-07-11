
import glob
import re

def replace_with_mx(x):
    print(x.group(0))
    libmex_path = x.group(1).split(';')[-1]
    libmx_path = libmex_path.replace("libmex.lib", "libmx.lib")
    replace_text = x.group(0) + ";" + libmx_path
    print(replace_text)
    return replace_text


def fix():

    files = glob.glob("build/mex_*.vcxproj")

    print("patching files: ", files)


    for f in files:

        print("patching: ", f)

        with open(f, 'r') as fh:
            file_string = fh.read()

        fixed_string = re.sub(r'<AdditionalDependencies>(.*?libmex.lib)', replace_with_mx, file_string, flags=re.DOTALL)

        with open(f, 'w') as fh:
            fh.write(fixed_string)







if __name__ == "__main__":
    fix()
