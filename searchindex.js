Search.setIndex({"docnames": ["floating_point", "index", "interpolation", "root_finding"], "filenames": ["floating_point.md", "index.md", "interpolation.md", "root_finding.md"], "titles": ["Floating Point Arithmetic", "Before you start", "Interpolation", "Root Finding"], "terms": {"In": [0, 2, 3], "thi": [0, 1, 2, 3], "chapter": [0, 3], "we": [0, 2, 3], "introduc": [0, 2], "some": [0, 2, 3], "basic": 0, "system": [0, 2], "modern": [0, 2], "discuss": [0, 2, 3], "ani": [0, 2, 3], "nonzero": [0, 2], "x": [0, 2, 3], "mathbb": [0, 2], "r": [0, 1, 2], "can": [0, 2, 3], "accur": 0, "repres": [0, 2], "an": [0, 1, 2, 3], "infinit": [0, 2], "sequenc": [0, 3], "digit": [0, 3], "understood": 0, "consequ": 0, "ration": 0, "ar": [0, 1, 2, 3], "interv": [0, 2, 3], "between": [0, 2, 3], "two": [0, 2, 3], "distinct": [0, 2], "i": [0, 1, 2, 3], "alwai": [0, 2, 3], "It": [0, 1, 2, 3], "fundament": [0, 1, 3], "properti": [0, 2], "therefor": [0, 2, 3], "binari": [0, 3], "write": [0, 2], "pm": [0, 2], "0": [0, 2, 3], "d_1": 0, "d_2": 0, "d_3": 0, "dot": [0, 2], "d_": 0, "t": [0, 2], "d_t": 0, "time": [0, 2], "e": [0, 2], "where": [0, 2, 3], "integ": 0, "expon": 0, "other": [0, 1, 2, 3], "d_i": 0, "The": [0, 1, 2, 3], "mantissa": 0, "note": [0, 3], "us": [0, 1, 2, 3], "differ": [0, 2, 3], "which": [0, 2, 3], "shift": 0, "left": [0, 2, 3], "one": [0, 2, 3], "bit": [0, 3], "There": [0, 2, 3], "essenti": [0, 3], "latter": [0, 3], "common": [0, 3], "practic": [0, 1, 3], "frac": [0, 2, 3], "cdot": [0, 2, 3], "order": [0, 2], "guarante": [0, 2, 3], "uniqu": [0, 2], "abov": [0, 2, 3], "need": [0, 2, 3], "further": 0, "assumpt": [0, 3], "exist": [0, 2, 3], "subset": 0, "n": [0, 1, 2, 3], "d_j": 0, "neq": [0, 2], "j": [0, 2, 3], "For": [0, 2, 3], "exampl": [0, 2, 3], "under": [0, 2], "111": 0, "take": [0, 2, 3], "gener": [0, 2, 3], "refer": 0, "set": [0, 2, 3], "finit": [0, 2], "length": 0, "consid": [0, 2, 3], "f": [0, 1, 2, 3], "e_": 0, "min": 0, "max": [0, 2], "mid": [0, 2], "le": [0, 2, 3], "cup": 0, "seen": [0, 2], "onli": [0, 2, 3], "smallest": 0, "posit": [0, 2], "element": [0, 2], "x_": [0, 2, 3], "largest": 0, "overlin": [0, 2], "call": [0, 2, 3], "normal": 0, "If": [0, 2, 3], "allow": [0, 3], "definit": [0, 3], "denorm": 0, "equidist": [0, 3], "distanc": 0, "proof": [0, 2, 3], "cap": 0, "d_1d_2": 0, "equidistantli": 0, "To": [0, 2], "understand": [0, 2], "approxim": [0, 1, 2, 3], "import": [0, 2], "maxim": 0, "rel": [0, 3], "respect": [0, 2], "closest": 0, "follow": [0, 1, 2, 3], "quantiti": 0, "max_": [0, 2], "min_": [0, 2], "z": [0, 2], "hold": [0, 2, 3], "mathrm": [0, 2], "u": [0, 1, 2, 3], "also": [0, 1, 2, 3], "unit": 0, "ha": [0, 2, 3], "version": [0, 3], "formal": 0, "appear": 0, "mostli": 0, "research": 0, "literatur": 0, "numer": [0, 1, 3], "packag": [0, 3], "lapack": 0, "program": [0, 1, 3], "languag": [0, 1, 3], "like": [0, 2], "python": [0, 1, 3], "matlab": [0, 1, 3], "c": [0, 1, 2, 3], "defin": [0, 2], "instead": [0, 3], "next": [0, 2, 3], "word": [0, 2], "correspond": [0, 2, 3], "strategi": 0, "former": [0, 3], "nearest": 0, "while": [0, 2, 3], "chop": 0, "necessari": [0, 2], "distinguish": [0, 2], "sinc": [0, 2, 3], "factor": 0, "without": [0, 2, 3], "loss": [0, 2, 3], "from": [0, 2, 3], "theorem": [0, 2, 3], "find": [0, 1, 2], "ast": [0, 3], "On": [0, 2, 3], "_": [0, 2, 3], "32": 0, "24": [0, 3], "125": 0, "128": 0, "quad": [0, 2, 3], "64": 0, "53": 0, "1021": 0, "1024": 0, "support": [0, 1, 3], "thei": 0, "often": 0, "singl": [0, 2], "doubl": 0, "slightli": [0, 2, 3], "instanc": [0, 2, 3], "underflow": 0, "all": [0, 2], "zero": [0, 2, 3], "float32": 0, "23": [0, 3], "sign": [0, 2, 3], "occupi": 0, "8": [0, 2, 3], "rang": 0, "126": 0, "127": 0, "store": [0, 2], "text": [0, 2, 3], "e_7": 0, "e_6": 0, "e_0": 0, "sum_": [0, 2, 3], "e_j": 0, "notic": [0, 2], "actual": [0, 2, 3], "have": [0, 2, 3], "256": 0, "valu": [0, 1, 2, 3], "reserv": 0, "special": [0, 2], "ones": [0, 3], "inf": 0, "infti": [0, 2, 3], "nan": 0, "Not": 0, "whether": 0, "lead": [0, 2, 3], "equal": [0, 2, 3], "quiet": 0, "signal": 0, "subnorm": 0, "occur": [0, 2, 3], "when": [0, 2, 3], "float64": 0, "similar": [0, 2, 3], "wai": [0, 2, 3], "11": [0, 2, 3], "52": 0, "textrm": 0, "fl": 0, "map": [0, 2], "Such": [0, 2], "written": [0, 2, 3], "out": 0, "explicitli": 0, "let": [0, 2, 3], "begin": [0, 2, 3], "case": [0, 2, 3], "end": [0, 2, 3], "clear": [0, 2, 3], "monoton": 0, "idempot": 0, "y": [0, 2], "rightarrow": 0, "trivial": [0, 2], "tild": [0, 2], "d": [0, 1, 2, 3], "_t": 0, "delta": [0, 2], "given": [0, 2, 3], "circ": 0, "outcom": 0, "straightforward": [0, 2], "realiz": 0, "box": 0, "divis": 0, "assum": [0, 2, 3], "corollari": [0, 2], "estim": [0, 2, 3], "align": [0, 2, 3], "close": [0, 2, 3], "magnitud": 0, "opposit": [0, 3], "signific": [0, 3], "deriv": [0, 2, 3], "complic": 0, "toward": [0, 3], "final": 0, "result": [0, 2, 3], "do": [0, 3], "techniqu": [0, 3], "fuse": 0, "multipli": 0, "add": [0, 2], "fma": 0, "here": [0, 2, 3], "quantifi": [0, 2, 3], "effect": [0, 2], "lemma": [0, 2], "a_1": [0, 2], "a_2": 0, "a_n": [0, 2], "a_k": [0, 2], "k": [0, 2, 3], "prod_": [0, 2], "b_n": [0, 2, 3], "quit": [0, 3], "easi": [0, 3], "induct": [0, 2], "b_1": [0, 2, 3], "suppos": [0, 2], "claim": 0, "m": [0, 1, 2, 3], "could": [0, 2], "see": [0, 2, 3], "b_m": 0, "a_": [0, 2], "b_": [0, 3], "impli": [0, 2, 3], "bound": [0, 2, 3], "naiv": [0, 2], "product": [0, 2], "p_n": 0, "x_j": [0, 2], "2n": [0, 2], "iter": 0, "p_k": 0, "x_1": [0, 2, 3], "p_": 0, "x_k": [0, 2], "ge": [0, 2, 3], "tau_k": 0, "th": [0, 2, 3], "step": [0, 3], "x_n": [0, 2, 3], "delta_n": 0, "delta_": [0, 2], "x_2": [0, 2, 3], "delta_j": 0, "tau_j": 0, "right": [0, 2, 3], "eta_n": 0, "s_n": 0, "s_k": [0, 2], "s_": 0, "carri": 0, "analysi": [0, 1, 2, 3], "_j": [0, 2], "denot": [0, 2, 3], "s_j": [0, 2], "l": [0, 2], "previou": [0, 2, 3], "provid": [0, 1, 2], "long": [0, 2, 3], "mapsto": [0, 2], "epsilon": [0, 3], "briefli": 0, "sum": 0, "reduct": 0, "A": [0, 1, 2, 3], "evalu": [0, 2, 3], "calcul": [0, 2, 3], "underbrac": [0, 2], "x_3": [0, 2], "x_4": [0, 2], "y_1": [0, 2], "y_2": 0, "y_3": 0, "y_i": 0, "obvious": 0, "y_": [0, 2], "t_n": 0, "prove": [0, 2], "nativ": [0, 3], "textup": 0, "h": [0, 2], "lceil": [0, 3], "log_2": [0, 2, 3], "rceil": [0, 3], "polynomi": 0, "p": [0, 2, 3], "a_0": [0, 2], "nest": 0, "form": [0, 3], "pleas": 0, "upper": [0, 3], "perimet": 0, "regular": 0, "polygon": 0, "inscrib": 0, "circumscrib": 0, "circl": 0, "diamet": 0, "start": [0, 2, 3], "hexagon": 0, "p_0": 0, "sqrt": [0, 2, 3], "equival": [0, 2, 3], "dfrac": 0, "each": [0, 2, 3], "side": [0, 3], "lim_": [0, 2, 3], "implement": [0, 3], "both": [0, 1, 2], "compar": 0, "exact": [0, 2, 3], "explain": [0, 3], "base": [0, 2, 3], "algorithm": [0, 2, 3], "mathcal": [0, 2], "o": 0, "b": [0, 2, 3], "keep": [0, 3], "track": 0, "input": [0, 3], "output": [0, 3], "get": [0, 3], "initi": [0, 1, 3], "y_j": [0, 2], "remov": 0, "perform": [0, 3], "restor": 0, "describ": [0, 3], "accuraci": 0, "test": 0, "randomli": 0, "sai": 0, "sim": [0, 2], "growth": 0, "expect": 0, "total": [0, 2, 3], "explan": 0, "your": 0, "you": 0, "higham": [0, 1], "1993": 0, "mari": 0, "2019": [0, 1], "muller": 0, "2006": 0, "hig93": 0, "nichola": 0, "siam": 0, "journal": 0, "scientif": [0, 2], "14": [0, 3], "783": 0, "799": 0, "hm19": 0, "theo": 0, "new": [0, 2, 3], "approach": [0, 3], "probabilist": 0, "41": 0, "a2815": 0, "a2835": 0, "mul06": 0, "jean": 0, "michel": 0, "elementari": 0, "function": [0, 2, 3], "springer": 0, "repositori": 1, "host": 1, "cours": 1, "materi": 1, "math": 1, "5630": 1, "6630": 1, "introduct": 1, "auburn": 1, "univers": 1, "textbook": 1, "first": [1, 2, 3], "method": [1, 2], "ascher": 1, "greif": 1, "2011": 1, "book": 1, "recommend": 1, "mathemat": [1, 3], "quarteroni": 1, "saaco": 1, "saleri": 1, "2007": 1, "theori": [1, 3], "lloyd": 1, "trefethen": 1, "ordinari": 1, "differenti": [1, 2, 3], "equat": [1, 3], "problem": [1, 2], "griffith": 1, "2010": 1, "concis": 1, "plato": 1, "2003": 1, "comput": [1, 2], "ieee": 1, "float": 1, "point": [1, 2, 3], "arithmet": 1, "overturn": 1, "2001": 1, "cover": 1, "topic": [1, 2], "root": [1, 2], "interpol": [1, 3], "integr": [1, 2], "solut": [1, 3], "design": [1, 2], "solid": 1, "foundat": 1, "applic": [1, 2], "requir": [1, 2, 3], "throughout": 1, "prerequisit": 1, "part": [1, 2], "linear": [1, 3], "2650": 1, "algebra": 1, "2660": 1, "default": 1, "class": 1, "script": 1, "julia": 1, "1d": 2, "solv": [2, 3], "type": [2, 3], "predefin": 2, "One": 2, "aid": 2, "cad": 2, "extens": 2, "manufactur": 2, "industri": 2, "speak": 2, "determin": 2, "paramet": 2, "access": 2, "pi_m": 2, "degre": 2, "seek": 2, "constraint": 2, "pi_n": 2, "y_k": 2, "1": 2, "node": 2, "resp": 2, "underdetermin": 2, "overdetermin": 2, "construct": 2, "a_j": 2, "formul": 2, "coeffici": 2, "pmatrix": 2, "x_0": [2, 3], "2": [2, 3], "vdot": 2, "ddot": 2, "y_0": 2, "y_n": 2, "matrix": 2, "v": 2, "vandermond": 2, "invert": 2, "Its": 2, "exercis": 2, "det": 2, "x_i": [2, 3], "g": [2, 3], "satisfi": [2, 3], "condit": 2, "most": [2, 3], "contradict": 2, "howev": [2, 3], "easier": 2, "somewhat": 2, "invers": 2, "l_0": 2, "l_1": 2, "l_n": 2, "l_j": 2, "jk": 2, "linearli": 2, "independ": 2, "basi": 2, "dimension": 2, "space": 2, "check": [2, 3], "preliminari": 2, "procedur": 2, "constant": [2, 3], "k_j": 2, "q": 2, "co": 2, "flop": 2, "advantag": [2, 3], "scheme": [2, 3], "anoth": [2, 3], "re": 2, "them": 2, "disadvantag": 2, "updat": [2, 3], "addit": 2, "cost": 2, "later": [2, 3], "overcom": 2, "issu": 2, "data": 2, "pair": [2, 3], "suffici": [2, 3], "smooth": 2, "possibl": [2, 3], "Then": [2, 3], "omega": 2, "xi": [2, 3], "roll": 2, "select": 2, "psi": 2, "chosen": [2, 3], "By": [2, 3], "least": [2, 3], "uniformli": 2, "number": [2, 3], "converg": 2, "interest": [2, 3], "think": 2, "about": 2, "convers": 2, "what": 2, "kind": 2, "vanish": 2, "tend": 2, "infin": 2, "10": [2, 3], "depend": 2, "size": [2, 3], "three": [2, 3], "term": [2, 3], "grow": 2, "rapidli": 2, "3": [2, 3], "larg": [2, 3], "so": [2, 3], "decai": 2, "show": [2, 3], "anymor": 2, "still": [2, 3], "true": [2, 3], "certain": [2, 3], "choic": [2, 3], "try": [2, 3], "better": 2, "difficult": 2, "worst": 2, "locat": [2, 3], "sub": 2, "sup_": 2, "4": [2, 3], "thu": [2, 3], "uniform": 2, "henc": 2, "4n": 2, "exponenti": 2, "work": 2, "awai": [2, 3], "origin": 2, "faster": [2, 3], "than": [2, 3], "would": [2, 3], "diverg": 2, "increas": 2, "famou": 2, "made": [2, 3], "carl": 2, "5": [2, 3], "shown": [2, 3], "around": 2, "6": [2, 3], "maximum": 2, "f_n": 2, "analyt": [2, 3], "domain": 2, "possibli": 2, "contain": 2, "hole": 2, "residu": 2, "simpl": [2, 3], "pole": 2, "pi": 2, "int_": 2, "partial": 2, "focus": 2, "studi": [2, 3], "behavior": 2, "equispac": 2, "over": [2, 3], "exp": 2, "int_a": 2, "log": 2, "sigma_n": 2, "contour": 2, "c_": [2, 3], "rho": [2, 3], "These": 2, "level": 2, "concentr": 2, "curv": 2, "midpoint": [2, 3], "enclos": 2, "insid": 2, "modulu": 2, "principl": 2, "must": 2, "attain": 2, "its": [2, 3], "boundari": 2, "small": 2, "situat": 2, "isol": 2, "z_k": 2, "rho_k": 2, "bigcup_": 2, "gamma_k": 2, "path": 2, "surround": 2, "ref": 2, "lem": 2, "ana": 2, "uni": 2, "con": 2, "whose": [2, 3], "limit": [2, 3], "goe": 2, "second": [2, 3], "summat": 2, "u_1": 2, "u_2": 2, "u_l": 2, "consist": [2, 3], "rho_1": 2, "rho_2": 2, "rho_m": 2, "decompos": 2, "group": 2, "substack": 2, "u_": 2, "As": 2, "doe": [2, 3], "whole": 2, "blow": 2, "up": 2, "why": [2, 3], "otherwis": [2, 3], "should": 2, "violat": 2, "intersect": 2, "real": 2, "line": [2, 3], "x_c": 2, "approx": [2, 3], "6334": 2, "onc": [2, 3], "pass": 2, "through": [2, 3], "prevent": 2, "section": [2, 3], "aim": 2, "minim": 2, "natur": 2, "question": 2, "restrict": 2, "our": 2, "simplic": [2, 3], "chang": [2, 3], "observ": [2, 3], "mean": [2, 3], "t_k": 2, "arcco": 2, "theta": [2, 3], "t_0": 2, "equiv": 2, "t_1": 2, "t_": 2, "extrema": 2, "statement": 2, "after": [2, 3], "replac": [2, 3], "variabl": [2, 3], "fourth": 2, "immedi": 2, "recurs": 2, "formula": [2, 3], "last": [2, 3], "more": 2, "importantli": 2, "optim": 2, "z_0": 2, "z_1": 2, "z_": 2, "due": [2, 3], "cancel": 2, "hand": [2, 3], "becaus": [2, 3], "share": 2, "same": [2, 3], "z_j": 2, "2j": 2, "now": [2, 3], "affin": 2, "phi": 2, "renewcommand": 2, "ep": 2, "varepsilon": 2, "much": 2, "smaller": [2, 3], "perturb": 2, "eps_j": 2, "_n": 2, "lambda_n": 2, "lebesgu": 2, "piecewis": 2, "lambda_": 2, "en": 2, "been": [2, 3], "paul": 2, "erd\u00f6": 2, "1964": 2, "faber": 2, "continu": [2, 3], "abl": 2, "almost": [2, 3], "sens": [2, 3], "extend": 2, "dynam": [2, 3], "scenario": 2, "alreadi": [2, 3], "found": 2, "f_k": 2, "how": [2, 3], "f_": 2, "automat": 2, "care": 2, "produc": [2, 3], "c_0": 2, "c_1": [2, 3], "c_2": [2, 3], "c_j": 2, "c_k": 2, "known": 2, "horner": 2, "c_3": 2, "innermost": 2, "c_n": [2, 3], "complex": 2, "3n": 2, "cheap": 2, "roughli": [2, 3], "5n": 2, "divid": 2, "squar": 2, "bracket": 2, "usual": [2, 3], "graph": 2, "help": 2, "relationship": 2, "among": 2, "searrow": 2, "main": 2, "power": 2, "g_n": 2, "7": [2, 3], "repeat": [2, 3], "moreov": 2, "taylor": [2, 3], "expans": [2, 3], "effici": [2, 3], "column": 2, "diagon": 2, "leadsto": 2, "color": 2, "red": 2, "green": 2, "cyan": 2, "blue": 2, "black": 2, "avail": [2, 3], "tupl": 2, "m_j": 2, "h_": 2, "pi_": 2, "idea": 2, "borrow": 2, "l_": 2, "x_l": 2, "obtain": [2, 3], "conclud": 2, "simplest": [2, 3], "diagram": 2, "fill": 2, "arrang": 2, "m_0": 2, "m_1": 2, "m_n": 2, "period": 2, "mani": [2, 3], "planar": 2, "parameter": 2, "suit": 2, "eventu": 2, "go": 2, "emph": 2, "kx": 2, "b_k": 2, "sin": 2, "said": [2, 3], "concept": 2, "valid": 2, "f_1": [2, 3], "f_2": [2, 3], "f_l": 2, "ident": 2, "rewrit": 2, "ik": 2, "gamma_0": 2, "gamma_": 2, "ib_k": 2, "substitut": 2, "ix": 2, "nx": 2, "exactli": 2, "nyquist": 2, "sampl": [2, 3], "direct": 2, "pari": 2, "simpli": 2, "creat": 2, "l_k": 2, "someth": 2, "rescal": 2, "remain": 2, "split": 2, "2x": [2, 3], "computation": 2, "reus": 2, "previous": 2, "barycentr": 2, "non": 2, "achiev": 2, "though": 2, "9": [2, 3], "ikx_j": 2, "im": 2, "gamma_m": 2, "emploi": 2, "reduc": [2, 3], "discret": 2, "texttt": [2, 3], "dft": 2, "vector": 2, "mathbf": 2, "fft": 2, "exploit": 2, "symmetri": [2, 3], "conquer": 2, "separ": 2, "even": [2, 3], "odd": 2, "half": [2, 3], "oper": 2, "2m": 2, "4m": 2, "includ": [2, 3], "multipl": 2, "2f": [2, 3], "4f": 2, "standard": [2, 3], "routin": 2, "softwar": 2, "focu": 2, "isx": 2, "remind": 2, "twice": [2, 3], "dx": 2, "absolut": [2, 3], "ly": 2, "region": 2, "contribut": [2, 3], "rest": [2, 3], "symmetr": 2, "inequ": [2, 3], "13": [2, 3], "h\u00f6lder": 2, "dunham": 2, "jackson": 2, "1913": 2, "improv": [2, 3], "low": 2, "pi_k": 2, "dimens": 2, "interfac": 2, "impos": 2, "2l": 2, "subinterv": 2, "explicit": 2, "s_1": 2, "combin": 2, "hat": 2, "theta_j": 2, "directli": 2, "well": [2, 3], "somewher": 2, "alpha": 2, "tau": 2, "scienc": 3, "engin": 3, "nonlinear": 3, "intermedi": 3, "fact": 3, "either": 3, "compact": 3, "mai": 3, "suffer": 3, "ineffici": 3, "instruct": 3, "risk": 3, "stack": 3, "overflow": 3, "varieti": 3, "tail": 3, "recus": 3, "compil": 3, "featur": 3, "execut": 3, "neither": 3, "nor": 3, "process": 3, "until": 3, "met": 3, "gain": 3, "desir": 3, "toler": 3, "veri": 3, "rough": 3, "taken": 3, "guess": 3, "sgn": 3, "regula": 3, "falsi": 3, "latin": 3, "account": 3, "li": 3, "connect": 3, "sometim": 3, "except": 3, "hard": 3, "never": 3, "decreas": 3, "illustr": 3, "fig": 3, "endpoint": 3, "fix": 3, "forc": 3, "weight": 3, "consecut": 3, "adjust": 3, "lambda": 3, "control": 3, "halv": 3, "reset": 3, "f_0": 3, "latest": 3, "return": 3, "At": 3, "modifi": 3, "make": 3, "closer": 3, "sever": 3, "termin": 3, "ftol": 3, "atol": 3, "rtol": 3, "aforement": 3, "cubic": 3, "5943130163548496": 3, "error": 3, "tabl": 3, "0000000000000000000000000": 3, "94e": 3, "01": 3, "5000000000000000000000000": 3, "43e": 3, "02": 3, "7500000000000000000000000": 3, "56e": 3, "6250000000000000000000000": 3, "07e": 3, "5625000000000000000000000": 3, "18e": 3, "5937500000000000000000000": 3, "63e": 3, "04": 3, "6093750000000000000000000": 3, "51e": 3, "6015625000000000000000000": 3, "25e": 3, "03": 3, "5976562500000000000000000": 3, "34e": 3, "5957031250000000000000000": 3, "39e": 3, "5947265625000000000000000": 3, "14e": 3, "12": 3, "5942382812500000000000000": 3, "47e": 3, "05": 3, "5944824218750000000000000": 3, "69e": 3, "5943603515625000000000000": 3, "73e": 3, "15": 3, "5942993164062500000000000": 3, "37e": 3, "16": 3, "5943298339843750000000000": 3, "68e": 3, "17": 3, "5943145751953125000000000": 3, "06": 3, "18": 3, "5943069458007812500000000": 3, "19": 3, "5943107604980468750000000": 3, "26e": 3, "20": 3, "5943126678466796875000000": 3, "49e": 3, "07": 3, "21": 3, "5943136215209960937500000": 3, "05e": 3, "22": 3, "5943131446838378906250000": 3, "28e": 3, "5943129062652587890625000": 3, "10e": 3, "5943130254745483398437500": 3, "12e": 3, "09": 3, "4444444444444446418174266": 3, "50e": 3, "5621621621621617492792211": 3, "22e": 3, "5876913365185605364615640": 3, "62e": 3, "5929610854818996301673906": 3, "35e": 3, "5940374914642010395482430": 3, "76e": 3, "5942568846837747997824408": 3, "61e": 3, "5943015817106331866170876": 3, "5943106870264029950590157": 3, "33e": 3, "5943125418534931370118102": 3, "75e": 3, "5943129196954899384763849": 3, "67e": 3, "08": 3, "6153846153846154187760931": 3, "11e": 3, "5847750865051901669744439": 3, "54e": 3, "5941951587569969106539247": 3, "5944267005726100450146987": 3, "5943130084597889606357057": 3, "90e": 3, "bigfloat": 3, "extract": 3, "precis": 3, "900": 3, "": 3, "5941951587569973547431346": 3, "5944267005726091568362790": 3, "5943130163543197674869134": 3, "29e": 3, "5943130163553775879847763": 3, "5943130163548486777358448": 3, "65e": 3, "25": 3, "13e": 3, "38": 3, "55e": 3, "75": 3, "70e": 3, "113": 3, "30e": 3, "225": 3, "few": 3, "identifi": 3, "robust": 3, "fast": 3, "rate": 3, "weaker": 3, "averag": 3, "quadrat": 3, "b_2": 3, "represent": 3, "b_i": 3, "b_j": 3, "argument": 3, "arriv": 3, "geometr": 3, "lower": 3, "move": 3, "larger": 3, "ell": 3, "slope": 3, "zeta": 3, "impl": 3, "unadjust": 3, "fall": 3, "below": 3, "reflect": 3, "becom": 3, "current": 3, "complet": 3, "cycl": 3, "quick": 3, "asymptot": 3, "epsilon_": 3, "simeq": 3, "appli": 3, "rule": 3, "epsilon_i": 3, "give": 3, "everi": 3, "analyz": 3, "variant": 3, "wa": 3, "discov": 3, "subroutin": 3, "ferranti": 3, "author": 3, "inform": 3, "dowel": 3, "jarratt": 3, "1972": 3, "discoveri": 3, "57": 3, "remark": 3, "ll": 3, "6f": 3, "full": 3, "7p": 3, "expand": 3, "retain": 3, "drop": 3, "correct": 3, "relat": 3, "although": 3, "gone": 3, "shorter": 3, "accord": 3, "prefer": 3, "finish": 3, "brief": 3, "avoid": 3, "swap": 3, "option": 3, "3rd": 3, "loop": 3, "elimin": 3, "illnoi": 3, "cannot": 3, "dj72": 3, "503": 3, "508": 3}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"float": 0, "point": 0, "arithmet": 0, "represent": 0, "real": 0, "number": 0, "what": 0, "doe": 0, "dens": 0, "mean": 0, "ieee754": 0, "standard": 0, "distribut": 0, "machin": 0, "precis": 0, "more": [0, 3], "about": [0, 3], "round": 0, "oper": 0, "cancel": 0, "error": [0, 2], "accumul": 0, "multipl": 0, "addit": 0, "exercis": [0, 3], "theoret": [0, 3], "part": [0, 3], "problem": [0, 3], "1": [0, 3], "2": 0, "3": 0, "horner": 0, "": [0, 2], "scheme": 0, "comput": [0, 3], "4": 0, "archimed": 0, "formula": 0, "pi": 0, "5": 0, "pairwis": 0, "summat": 0, "6": 0, "kahan": 0, "compens": 0, "7": 0, "extend": [0, 3], "read": [0, 3], "befor": 1, "you": 1, "start": 1, "interpol": 2, "polynomi": 2, "lagrang": 2, "rung": 2, "phenomenon": 2, "remaind": 2, "theori": 2, "chebyshev": 2, "stabil": 2, "newton": [2, 3], "form": 2, "hermit": 2, "trigonometr": 2, "fourier": 2, "seri": 2, "fast": 2, "transform": 2, "spline": 2, "linear": 2, "root": 3, "find": 3, "bracket": 3, "method": 3, "bisect": 3, "iter": 3, "v": 3, "recurs": 3, "fals": 3, "posit": 3, "illinoi": 3, "stop": 3, "criteria": 3, "order": 3, "converg": 3, "big": 3, "o": 3, "notat": 3, "select": 3, "decai": 3, "factor": 3, "pegasu": 3, "anderson": 3, "bjorck": 3, "raphson": 3, "secant": 3, "applic": 3, "optim": 3}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 60}, "alltitles": {"Floating Point Arithmetic": [[0, "floating-point-arithmetic"]], "Representation of Real Numbers": [[0, "representation-of-real-numbers"]], "What does \u201cdense\u201d mean?": [[0, null]], "IEEE754 Standard": [[0, null]], "Floating Point Numbers": [[0, "floating-point-numbers"]], " (Distribution of Floating Numbers)": [[0, "THM-Di-Fl-Nu"]], " (Machine Precision)": [[0, "THM-Ma-Pr"]], "": [[0, null], [0, "THM-REL-ERR"], [0, "COR-REL-ERR"], [0, "LEM-ACC"], [0, "definition-6"], [2, "THM-INTER-UNIQ"], [2, "DEF-LA-PO"], [2, "THM-UNIQ-LAG-INTER"], [2, "remark-3"], [2, "THM-INTERP-ERROR"], [2, "corollary-5"], [2, "example-6"], [2, "Lem-2-Ome-Lim"], [2, "Lem-2-Ana-Uni-Con"], [2, "example-9"], [2, "definition-10"], [2, "theorem-11"], [2, "theorem-12"], [2, "definition-13"], [2, "COR-CHEBY"], [2, "remark-15"], [2, "definition-16"], [2, "theorem-17"], [2, "remark-18"], [2, "remark-19"], [2, "remark-20"], [2, "remark-21"], [2, "remark-22"], [2, "definition-23"], [2, "lemma-24"], [2, "remark-25"], [2, "corollary-26"], [2, "remark-27"], [2, "THM-TRIG"], [2, "theorem-29"], [2, "definition-30"], [3, null], [3, "example-1"], [3, "rmk:bracket-methods"], [3, "dfn-order-of-convergence"], [3, "thm-bisection-convergence"], [3, "thm-illinois-convergence"], [3, "rmk-illinois-convergence"], [3, "rmk-pegasus-asymptotic"], [3, null], [3, "thm-anderson-bjorck-convergence"]], "More About IEEE754 Standard": [[0, null]], "Rounding": [[0, "rounding"]], "Arithmetic Operations": [[0, "arithmetic-operations"]], " (Cancellation Error)": [[0, "remark-4"]], "Error Accumulation: Multiplication": [[0, "error-accumulation-multiplication"]], "Error Accumulation: Addition": [[0, "error-accumulation-addition"]], "Exercises": [[0, "exercises"], [3, "exercises"]], "Theoretical Part": [[0, "theoretical-part"], [3, "theoretical-part"]], "Problem 1": [[0, null], [3, null]], "Problem 2": [[0, null]], "Problem 3 (Horner\u2019s scheme)": [[0, null]], "Computational Part": [[0, "computational-part"], [3, "computational-part"]], "Problem 4 (Archimedes\u2019 formula for \\pi)": [[0, null]], "Problem 5 (Pairwise summation)": [[0, null]], "Problem 6 (Kahan compensated summation)": [[0, null]], " (Kahan Compensated Summation)": [[0, "AL-KA-CO-SU"]], "Problem 7": [[0, null]], "Extended Reading": [[0, "extended-reading"], [3, "extended-reading"]], "Before you start": [[1, "before-you-start"]], "Interpolation": [[2, "interpolation"]], "Polynomial Interpolation": [[2, "polynomial-interpolation"]], "Lagrange Polynomial": [[2, "lagrange-polynomial"]], "Interpolation Error": [[2, "interpolation-error"]], "Runge\u2019s Phenomenon": [[2, "runge-s-phenomenon"]], "Interpolation Remainder Theory": [[2, "interpolation-remainder-theory"]], "Chebyshev Interpolation": [[2, "chebyshev-interpolation"]], "Stability of Polynomial Interpolation": [[2, "stability-of-polynomial-interpolation"]], "Newton Form": [[2, "newton-form"]], "Hermite Polynomial Interpolation": [[2, "hermite-polynomial-interpolation"]], "Trigonometric Interpolation": [[2, "trigonometric-interpolation"]], "Fourier Series": [[2, "fourier-series"]], "Fast Fourier Transform": [[2, "fast-fourier-transform"]], "Interpolation Error of Trigonometric Polynomial": [[2, "interpolation-error-of-trigonometric-polynomial"]], "Spline Interpolation": [[2, "spline-interpolation"]], "Linear Spline": [[2, "linear-spline"]], "Root Finding": [[3, "root-finding"]], "Bracket Methods": [[3, "bracket-methods"]], "Bisection Method": [[3, "bisection-method"]], "Iterative vs Recursive": [[3, null]], "False Position Method": [[3, "false-position-method"]], " (Illinois Method)": [[3, "AL-ILLINOIS"]], "Stopping Criteria": [[3, null]], "Order of Convergence": [[3, "order-of-convergence"]], "More About Convergence": [[3, null]], "Big O Notation": [[3, null]], "Selection of Decay Factor": [[3, null]], " (Pegasus Method)": [[3, "thm-pegasus-convergence"]], " (Anderson-Bjorck Method)": [[3, "AL-ANDERSON-BJORCK"]], "Iterative Methods": [[3, "iterative-methods"]], "Newton-Raphson Method": [[3, "newton-raphson-method"]], "Secant Method": [[3, "secant-method"]], "Applications in Optimization": [[3, "applications-in-optimization"]]}, "indexentries": {}})