// Version dropdown, archived-version banner and version history for the AF
// formulations page. Shared by the live page and its frozen snapshots in
// versions/<id>/; everything comes from versions.json.
//   <script src="versions.js" data-version="latest" data-root="./"></script>        live page
//   <script src="../../versions.js" data-version="d3" data-root="../../"></script>  snapshot
(() => {
  const me = document.currentScript;
  const root = me.dataset.root || "./";
  const want = me.dataset.version || "latest";

  const css = `
.vbar { display: flex; flex-wrap: wrap; align-items: center; gap: 6px 10px; margin: -4px 0 14px; font-size: 0.85rem; color: var(--ink-2, #52514e); }
.vbar label { color: var(--muted, #898781); }
.vbar select { font: inherit; color: var(--ink, #0b0b0b); background: var(--surface, #fcfcfb); border: 1px solid var(--border, rgba(11,11,11,.1)); border-radius: 999px; padding: 3px 8px; max-width: 100%; min-width: 0; flex: 0 1 auto; }
.vbar a { color: var(--f3, #2a78d6); }
.vbanner { background: var(--warn-bg, #fff4d6); color: var(--warn-ink, #5a3d00); border-left: 4px solid var(--warn-border, #eda100); padding: 10px 14px; border-radius: 6px; margin: 12px 0 16px; font-size: 0.95rem; }
.vbanner a { color: inherit; font-weight: 600; }
.vhist { list-style: none; padding: 0; margin: 0; }
.vhist > li { border-top: 1px solid var(--grid, #e1e0d9); padding: 12px 0; }
.vhist > li:first-child { border-top: 0; }
.vhist .vh-head { display: flex; flex-wrap: wrap; gap: 4px 10px; align-items: baseline; }
.vhist .vh-meta { font-size: 0.8rem; color: var(--muted, #898781); }
.vhist .vh-flag { font-size: 0.75rem; border: 1px solid var(--warn-border, #eda100); color: var(--warn-ink, #5a3d00); background: var(--warn-bg, #fff4d6); border-radius: 999px; padding: 1px 8px; }
.vhist ul { margin: 6px 0 0; padding-left: 1.2em; }
.vhist li li { margin: 3px 0; }
`;
  const st = document.createElement("style"); st.textContent = css; document.head.appendChild(st);

  const fmtDate = d => new Date(d + "T12:00:00").toLocaleDateString("en-GB", { day: "numeric", month: "short", year: "numeric" });
  const href = (v, latest) => v.id === latest ? root : root + "versions/" + v.id + "/";
  const h = (tag, props = {}, kids = []) => { const e = document.createElement(tag); Object.assign(e, props); kids.forEach(k => e.append(k)); return e; };
  const anchor = () => document.querySelector("main h1") || document.querySelector("h1");

  fetch(root + "versions.json", { cache: "no-cache" }).then(r => r.json()).then(V => {
    const latest = V.latest, list = V.versions || [];
    const id = want === "latest" ? latest : want;
    const cur = list.find(v => v.id === id) || { id, date: "", title: id };
    const archived = id !== latest;

    // dropdown
    const sel = h("select", { id: "vsel", ariaLabel: "Page version" });
    list.forEach(v => sel.append(h("option", { value: v.id, selected: v.id === id,
      textContent: `${v.id === latest ? "Latest · " : ""}${v.id} · ${fmtDate(v.date)} · ${v.title.replace(/^Draft \d+: (.)/, (_, c) => c.toUpperCase())}` })));
    sel.addEventListener("change", () => { const v = list.find(x => x.id === sel.value); if (v) location.href = href(v, latest); });
    const bar = h("div", { className: "vbar" }, [h("label", { htmlFor: "vsel", textContent: "Version" }), sel,
      h("a", { href: root + "#history", textContent: "What changed between versions" })]);
    const a = anchor(); if (a) a.after(bar); else document.body.prepend(bar);

    // archived snapshot: banner; comments closed (CSS in the snapshot also hides them)
    if (archived) {
      const L = list.find(v => v.id === latest);
      bar.after(h("div", { className: "vbanner", role: "note" }, [
        h("strong", { textContent: `You are reading an earlier version: ${cur.title} (${fmtDate(cur.date)}). ` }),
        "Its numbers may have been revised since. Comments are closed on earlier versions. ",
        h("a", { href: root, textContent: `Read the latest version (${latest}${L ? ", " + fmtDate(L.date) : ""}) →` })]));
      document.querySelectorAll(".cbtn, #cdrawer").forEach(e => e.remove());
    }

    // version history (live page only)
    const box = document.getElementById("vhistory");
    if (box && archived) box.replaceChildren(h("p", {}, ["The current version history is on the ", h("a", { href: root + "#history", textContent: "latest page" }), "."]));
    if (box && !archived) {
      const ol = h("ol", { className: "vhist" });
      list.forEach(v => {
        const head = h("div", { className: "vh-head" }, [
          h("b", {}, [v.id === latest ? `${v.id} · ${v.title} (this page)` : h("a", { href: href(v, latest), textContent: `${v.id} · ${v.title}` })]),
          h("span", { className: "vh-meta", textContent: `${fmtDate(v.date)} · ${v.ref ? "tag " + v.ref : "commit " + v.commit}` })]);
        if (v.results_changed) head.append(h("span", { className: "vh-flag", textContent: "results changed" }));
        ol.append(h("li", {}, [head, h("p", { textContent: v.summary, style: "margin:4px 0 0" }),
          h("ul", {}, v.changes.map(c => h("li", { textContent: c })))]));
      });
      box.replaceChildren(ol);
    }
  }).catch(() => {
    const a = anchor(); if (a) a.after(h("div", { className: "vbar", textContent: "Version " + want }));
  });
})();
