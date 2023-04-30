# ###### Art divrtsion
# Some yout generative produced while figuring out the mappoing from (many) colorspaces to RGBW.
# Speficially for LEDs.  Was not trivial, and hinged on a friend ultimately working it out in an espteric soace, HSL.
#
#  github: https://github.com/iamh2o/rgbw_colorspace_converter


localrules:
    art,
    arta,


rule arta:  # TARGET: A visual derrive.
    params:
        ainfo="Hi-  the following  script will be run run_color_module_RGB_HSV_HEX_demo.py , it wil roll for a while before haltine, and usually will drpop a pndg and html file of o=your visuals in PWD. To exit you'll need to probably hammer on crrl-c...... and 50/50 odss when your screen returns the cursor will be  invisible.   No fear, hih enter a few times, then carefully type 'clear' -enter- then 'reset' enter and you'll be back to normal.   You can pass arguments to this script via the command line --config art_args=' -n -b___-//=. -u 40 -g -y -z -j' .  If no args are supplied, the bais mode willl run with a random seed    ---- enjoy ",
        art_args="-j" if "art_args" not in config else config["art_args"],
    output:
        "logs/made.art",
    shell:
        """
        if [[ "{params.art_args}" == "" ]]; then   \
            {params.ainfo};
            clear;
            reset;
            timeout 10  run_color_module_RGB_HSV_HEX_demo.py {params.art_args} ;  \
        else   \
            run_color_module_RGB_HSV_HEX_demo.py {params.art_args} ;   \
        fi; \
         \
        clear; \
        reset; \
        touch {output}; \
        """


rule art:  # TARGET: A visual derrive.
    input:
        "logs/made.art",
    shell:
        "rm {input};"
