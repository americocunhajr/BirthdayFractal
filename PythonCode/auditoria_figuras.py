"""
Auditoria das figuras do artigo: para cada figura de fractal, recalcula exatamente
com os parametros usados na renderizacao (grade 3200x3200, ate 300 iteracoes) e
conta os pontos de cada categoria:
  0 classificado (convergiu para uma raiz pelo criterio numerico);
  1 valor nao finito em ponto flutuante (passo ou iterado inf/nan: zero exato de p',
    forma 0/0 numa raiz multipla, transbordamento etc.; registro numerico, sem
    significado matematico proprio);
  2 nao classificado apos 300 iteracoes;
  3 parou (passo pequeno) longe de toda raiz.
Para a Figura 10 (z^3-2z+2), verifica tambem se os pontos nao classificados
terminaram perto do ciclo {0,1}.

Uso:  python auditoria_figuras.py            (todas as figuras, em sequencia)
      python auditoria_figuras.py i n        (parte i de n, para rodar em paralelo)
      python auditoria_figuras.py juntar n   (junta as n partes)
Saida: auditoria_figuras.json (ou auditoria_figuras_parte<i>.json).
"""
import sys
import json
import numpy as np

from newton_arte import calcular, janela_padrao
from newton_aniversario import coeficientes

DATAS_FIG = ["25/12/1642", "04/01/1643", "30/04/1777", "15/04/1707", "16/08/1821",
             "29/04/1854", "28/02/1878", "03/02/1893", "20/11/1924", "05/07/1687",
             "23/09/1846", "30/11/1889", "29/05/1919", "29/12/1972", "16/05/1980",
             "08/08/1900", "30/08/1984", "03/07/1986", "27/04/1952", "15/08/1950"]

casos = []
for d in DATAS_FIG:
    casos.append((f"newton_{d.replace('/', '-')}", coeficientes(d), {}))
casos.append(("cayley_z2-1", np.array([1., 0, -1]), {}))
casos.append(("cayley_z3-1", np.array([1., 0, 0, -1]), {}))
cf = np.array([1., 0, -2, 2])
cenf, Rf = janela_padrao(cf)
casos.append(("falha_z3-2z+2", cf, dict(centro=0.5 + 0j, R=1.0, Rref=Rf)))
cg = np.array(coeficientes("30/04/1777"), float)
cen, R = janela_padrao(cg)
casos.append(("zoom_2", cg, dict(centro=complex(-0.45128, 0.0), R=R / 6, Rref=R)))
casos.append(("zoom_3", cg, dict(centro=complex(-0.6158, 0.0919), R=R / 48, Rref=R)))



def auditar(lista):
    res = {}
    for nome, c, kw in lista:
        bacia, nu, raizes, centro, Rj, info = calcular(np.asarray(c, float), n=1600, ss=2,
                                                        com_status=True, **kw)
        st = info["status"]
        cont = {str(k): int(np.sum(st == k)) for k in range(4)}
        item = {"pontos": int(st.size), "contagem": cont,
                "centro": [float(np.real(centro)), float(np.imag(centro))],
                "meia_largura": float(Rj)}
        if nome == "falha_z3-2z+2":
            Zf = info["Zfinal"][st == 2]
            perto = np.minimum(np.abs(Zf - 0), np.abs(Zf - 1))
            item["nao_classificados_perto_do_ciclo_1e-6"] = int(np.sum(perto < 1e-6))
        res[nome] = item
        print(nome, cont, item.get("nao_classificados_perto_do_ciclo_1e-6", ""), flush=True)
    return res


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        json.dump(auditar(casos), open("auditoria_figuras.json", "w"), indent=1)
    elif a[0] == "juntar":
        res = {}
        for i in range(int(a[1])):
            res.update(json.load(open(f"auditoria_figuras_parte{i}.json")))
        ordem = [nome for nome, _, _ in casos]
        assert sorted(res) == sorted(ordem), "partes incompletas"
        json.dump({k: res[k] for k in ordem}, open("auditoria_figuras.json", "w"), indent=1)
        print("auditoria_figuras.json com", len(res), "figuras")
    else:
        i, n = int(a[0]), int(a[1])
        json.dump(auditar(casos[i::n]), open(f"auditoria_figuras_parte{i}.json", "w"), indent=1)
