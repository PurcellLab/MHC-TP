/* Adds a "Copy page as text" button to every docs page. Copies the rendered
   text of the main content to the clipboard. Re-runs on Material instant
   navigation. */
(function () {
  function addCopyButton() {
    var content = document.querySelector(".md-content__inner");
    if (!content || content.querySelector(".copy-page-btn")) return;

    var btn = document.createElement("button");
    btn.type = "button";
    btn.className = "copy-page-btn";
    btn.innerHTML = '<span class="cpb-label">Copy page as text</span>';

    btn.addEventListener("click", function () {
      var clone = content.cloneNode(true);
      clone.querySelectorAll(".copy-page-btn").forEach(function (b) { b.remove(); });
      var text = clone.innerText.replace(/\n{3,}/g, "\n\n").trim();

      function done() {
        var label = btn.querySelector(".cpb-label");
        var prev = label.textContent;
        label.textContent = "Copied!";
        btn.classList.add("copied");
        setTimeout(function () {
          label.textContent = prev;
          btn.classList.remove("copied");
        }, 1500);
      }

      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(done).catch(fallback);
      } else {
        fallback();
      }

      function fallback() {
        var ta = document.createElement("textarea");
        ta.value = text;
        ta.style.position = "fixed";
        ta.style.opacity = "0";
        document.body.appendChild(ta);
        ta.focus();
        ta.select();
        try { document.execCommand("copy"); done(); } catch (e) { /* ignore */ }
        ta.remove();
      }
    });

    content.insertBefore(btn, content.firstChild);
  }

  if (typeof window.document$ !== "undefined" && window.document$.subscribe) {
    window.document$.subscribe(addCopyButton); // Material instant navigation
  } else {
    document.addEventListener("DOMContentLoaded", addCopyButton);
  }
})();
